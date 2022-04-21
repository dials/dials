# LIBTBX_SET_DISPATCHER_NAME dials.refine_bravais_settings
# LIBTBX_SET_DISPATCHER_NAME dials.rbs

"""Refinement of Bravais settings consistent with the primitive unit cell.

This program takes as input the output of dials.index, i.e. indexed.expt
and indexed.refl files. Full refinement of the crystal and experimental
geometry parameters will be performed (by default) in all Bravais settings
that are consistent with the input primitive unit cell. A table is printed
containing various information for each potential Bravais setting, including
the metric fit (a measure of the deviation from the triclinic cell),
the root-mean-square-deviations (rmsd), in mm, between the observed and
predicted spot centroids, the refined unit cell parameters in each Bravais
setting, and the change of basis operator to transform from the triclinic cell
to each Bravais setting.

The program also generates a .expt file for each Bravais setting, e.g.
bravais_setting_1.expt, which is equivalent to the input indexed.expt, but
with the crystal model refined in the chosen Bravais setting. These
bravais_setting_*.expt files are suitable as input to dials.refine or
dials.integrate, although the indexed.refl file will need to be re-indexed
using dials.reindex if the change of basis operator (cb_op) for the chosen
Bravais setting is not the identity operator (a,b,c).

Examples::

  dials.refine_bravais_settings indexed.expt indexed.refl

  dials.refine_bravais_settings indexed.expt indexed.refl nproc=4
"""


from __future__ import annotations

import collections
import json
import logging
import os
import sys

import libtbx.phil
from cctbx import sgtbx
from cctbx.sgtbx import bravais_types
from dxtbx.model import ExperimentList

import dials.util
from dials.algorithms.indexing.bravais_settings import (
    refined_settings_from_refined_triclinic,
)
from dials.array_family import flex
from dials.util import log
from dials.util.options import ArgumentParser, reflections_and_experiments_from_files
from dials.util.version import dials_version

logger = logging.getLogger("dials.command_line.refine_bravais_settings")

phil_scope = libtbx.phil.parse(
    """
include scope dials.algorithms.indexing.bravais_settings.phil_scope

crystal_id = None
  .type = int(value_min=0)

output {
  directory = "."
    .type = path
  log = dials.refine_bravais_settings.log
    .type = path
  prefix = None
    .type = str
}

""",
    process_includes=True,
)


def bravais_lattice_to_space_groups(chiral_only=True):
    bravais_lattice_to_sg = collections.defaultdict(list)
    for sgn in range(230):
        sg = sgtbx.space_group_info(number=sgn + 1).group()
        if (not chiral_only) or (sg.is_chiral()):
            bravais_lattice = bravais_types.bravais_lattice(group=sg)
            bravais_lattice_to_sg[str(bravais_lattice)].append(sg)
    bravais_lattice_to_sg["mI"] = [sgtbx.space_group_info("I2").group()]
    return bravais_lattice_to_sg


def bravais_lattice_to_space_group_table(bravais_settings=None, chiral_only=True):
    bravais_lattice_to_sg = bravais_lattice_to_space_groups(chiral_only=chiral_only)
    logger.info("Chiral space groups corresponding to each Bravais lattice:")
    for bravais_lattice, space_groups in bravais_lattice_to_sg.items():
        if bravais_settings is not None and bravais_lattice not in bravais_settings:
            continue
        logger.info(
            ": ".join(
                [
                    bravais_lattice,
                    " ".join([short_space_group_name(sg) for sg in space_groups]),
                ]
            )
        )


def short_space_group_name(space_group):
    sgt = space_group.type()
    symbol = sgt.lookup_symbol()
    if sgt.number() > 1:
        symbol = symbol.replace(" 1", "")
    return symbol.replace(" ", "")


def eliminate_sys_absent(experiments, reflections):
    if experiments[0].crystal.get_space_group().n_ltr() > 1:
        effective_group = (
            experiments[0]
            .crystal.get_space_group()
            .build_derived_reflection_intensity_group(anomalous_flag=True)
        )
        sys_absent_flags = effective_group.is_sys_absent(reflections["miller_index"])
        reflections = reflections.select(~sys_absent_flags)
    return reflections


def map_to_primitive(experiments, reflections):
    """Map experiments and reflections to primitive setting."""
    cb_op_to_primitive = (
        experiments[0]
        .crystal.get_space_group()
        .info()
        .change_of_basis_op_to_primitive_setting()
    )
    experiments[0].crystal.update(
        experiments[0].crystal.change_basis(cb_op_to_primitive)
    )
    reflections["miller_index"] = cb_op_to_primitive.apply(reflections["miller_index"])


def select_datasets_on_crystal_id(experiments, reflections, crystal_id):
    """Select experiments and reflections with the given crystal id"""
    assert crystal_id < len(experiments.crystals())
    experiment_ids = experiments.where(crystal=experiments.crystals()[crystal_id])

    experiments = ExperimentList([experiments[i] for i in experiment_ids])
    refl_selections = [reflections["id"] == i for i in experiment_ids]
    reflections["id"] = flex.int(len(reflections), -1)
    for i, sel in enumerate(refl_selections):
        reflections["id"].set_selected(sel, i)
    reflections = reflections.select(reflections["id"] > -1)
    return experiments, reflections


@dials.util.show_mail_handle_errors()
def run(args=None):
    usage = "dials.refine_bravais_settings indexed.expt indexed.refl [options]"

    parser = ArgumentParser(
        usage=usage,
        phil=phil_scope,
        read_experiments=True,
        read_reflections=True,
        check_format=False,
        epilog=__doc__,
    )

    params, options = parser.parse_args(args=args, show_diff_phil=False)

    # Configure the logging
    log.config(verbosity=options.verbose, logfile=params.output.log)

    logger.info(dials_version())

    # Log the diff phil
    diff_phil = parser.diff_phil.as_str()
    if diff_phil != "":
        logger.info("The following parameters have been modified:\n")
        logger.info(diff_phil)

    reflections, experiments = reflections_and_experiments_from_files(
        params.input.reflections, params.input.experiments
    )
    if len(reflections) == 0 or len(experiments) == 0:
        parser.print_help()
        return

    assert len(reflections) == 1
    reflections = reflections[0]

    # Reduce what we pickle and send to workers by removing unused data
    if "shoebox" in reflections:
        del reflections["shoebox"]

    if len(experiments) == 0:
        parser.print_help()
        return
    if len(experiments.crystals()) > 1:
        if params.crystal_id is not None:
            experiments, reflections = select_datasets_on_crystal_id(
                experiments, reflections, params.crystal_id
            )
        else:
            sys.exit(
                "Only one crystal can be processed at a time: set crystal_id to choose "
                "experiment."
            )

    reflections = eliminate_sys_absent(experiments, reflections)
    map_to_primitive(experiments, reflections)

    refined_settings = refined_settings_from_refined_triclinic(
        experiments, reflections, params
    )
    possible_bravais_settings = {solution["bravais"] for solution in refined_settings}
    bravais_lattice_to_space_group_table(possible_bravais_settings)
    logger.info(refined_settings)

    prefix = params.output.prefix
    if prefix is None:
        prefix = ""
    summary_file = f"{prefix}bravais_summary.json"
    logger.info("Saving summary as %s", summary_file)
    with open(os.path.join(params.output.directory, summary_file), "w") as fh:
        json.dump(refined_settings.as_dict(), fh)

    for subgroup in refined_settings:
        expts = subgroup.refined_experiments
        soln = int(subgroup.setting_number)
        bs_json = "%sbravais_setting_%i.expt" % (prefix, soln)
        logger.info("Saving solution %i as %s", soln, bs_json)
        expts.as_file(os.path.join(params.output.directory, bs_json))


if __name__ == "__main__":
    run()

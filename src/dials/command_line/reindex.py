# DIALS_ENABLE_COMMAND_LINE_COMPLETION


from __future__ import annotations

import logging
import os
import sys

import numpy as np

import iotbx.phil
from cctbx import sgtbx
from libtbx import Auto

import dials.util
from dials.algorithms.indexing.assign_indices import AssignIndicesGlobal
from dials.algorithms.scaling.scaling_library import determine_best_unit_cell
from dials.array_family import flex
from dials.util import Sorry, log
from dials.util.filter_reflections import filtered_arrays_from_experiments_reflections
from dials.util.options import ArgumentParser, reflections_and_experiments_from_files
from dials.util.reference import intensities_from_reference_file
from dials.util.reindex import (
    change_of_basis_op_against_reference,
    derive_change_of_basis_op,
    reindex_experiments,
    reindex_reflections,
)
from dials.util.version import dials_version

logger = logging.getLogger("dials.command_line.reindex")

help_message = """

This program can be used to re-index an indexed.expt and/or indexed.refl
file from one setting to another. The change of basis operator can be
provided in h,k,l, or a,b,c or x,y,z conventions. By default the change of
basis operator will also be applied to the space group in the indexed.expt
file, however, optionally, a space group (including setting) to be applied
AFTER applying the change of basis operator can be provided.
Alternatively, to reindex an integrated dataset in the case of indexing ambiguity,
a reference dataset (models.expt and reflection.refl) in the same space
group can be specified. In this case, any potential twin operators are tested,
and the dataset is reindexed to the setting that gives the highest correlation
with the reference dataset.

Examples::

  dials.reindex indexed.expt change_of_basis_op=b+c,a+c,a+b

  dials.reindex indexed.refl change_of_basis_op=-b,a+b+2*c,-a

  dials.reindex indexed.expt indexed.refl change_of_basis_op=l,h,k

  dials.reindex indexed.expt indexed.refl reference.experiments=reference.expt
    reference.reflections=reference.refl
"""

phil_scope = iotbx.phil.parse(
    """
change_of_basis_op = a,b,c
  .type = str
hkl_offset = None
  .type = ints(size=3)
space_group = None
  .type = space_group
  .help = "The space group to be applied AFTER applying the change of basis "
           "operator."
reference {
  experiments = None
    .type = path
    .help = "Reference experiment for determination of change of basis operator."
  reflections = None
    .type = path
    .help = "Reference reflections to allow reindexing to consistent index between datasets."
  file = None
    .type = path
    .help = "A file containing a reference set of intensities e.g. MTZ/cif, or a"
            "file from which a reference set of intensities can be calculated"
            "e.g. .pdb or .cif . The space group of the reference file will"
            "be used and if an indexing ambiguity is present, the input"
            "data will be reindexed to be consistent with the indexing mode of"
            "this reference file."
    .expert_level = 2
  include scope dials.util.reference.reference_phil_str
}
output {
  experiments = reindexed.expt
    .type = str
    .help = "The filename for reindexed experimental models"

  reflections = reindexed.refl
    .type = str
    .help = "The filename for reindexed reflections"
  log = dials.reindex.log
    .type = path
}
""",
    process_includes=True,
)


@dials.util.show_mail_handle_errors()
def run(args=None):
    usage = "dials.reindex [options] indexed.expt indexed.refl"

    parser = ArgumentParser(
        usage=usage,
        phil=phil_scope,
        read_reflections=True,
        read_experiments=True,
        check_format=False,
        epilog=help_message,
    )

    params, options = parser.parse_args(args, show_diff_phil=False)

    log.config(verbosity=options.verbose, logfile=params.output.log)

    logger.info(dials_version())

    diff_phil = parser.diff_phil.as_str()
    if diff_phil != "":
        logger.info("The following parameters have been modified:\n")
        logger.info(diff_phil)

    reflections, experiments = reflections_and_experiments_from_files(
        params.input.reflections, params.input.experiments
    )
    if len(experiments) == 0 and len(reflections) == 0:
        parser.print_help()
        return
    if params.change_of_basis_op is None:
        raise Sorry("Please provide a change_of_basis_op.")

    reference_crystal = None
    if params.reference.experiments is not None:
        from dxtbx.serialize import load

        reference_experiments = load.experiment_list(
            params.reference.experiments, check_format=False
        )
        if len(reference_experiments.crystals()) == 1:
            reference_crystal = reference_experiments.crystals()[0]
        else:
            # first check sg all same
            sgs = [
                expt.crystal.get_space_group().type().number() for expt in experiments
            ]
            if len(set(sgs)) > 1:
                raise Sorry(
                    """The reference experiments have different space groups:
                    space group numbers found: %s
                    Please reanalyse the data so that space groups are consistent,
                    (consider using dials.reindex, dials.symmetry or dials.cosym)"""
                    % ", ".join(map(str, set(sgs)))
                )

            reference_crystal = reference_experiments.crystals()[0]
            reference_crystal.set_unit_cell(
                determine_best_unit_cell(reference_experiments)
            )

    if params.reference.reflections is not None:
        # First check that we have everything as expected for the reference reindexing
        if params.reference.experiments is None:
            raise Sorry(
                """For reindexing against a reference dataset, a reference
experiments file must also be specified with the option: reference.experiments= """
            )
        if not os.path.exists(params.reference.reflections):
            raise Sorry("Could not locate reference dataset reflection file")

        reference_reflections = flex.reflection_table().from_file(
            params.reference.reflections
        )
        if (
            reference_reflections.get_flags(
                reference_reflections.flags.integrated_sum
            ).count(True)
            == 0
        ):
            assert "intensity.sum.value" in reference_reflections, (
                "No 'intensity.sum.value in reference reflections"
            )
            reference_reflections.set_flags(
                flex.bool(reference_reflections.size(), True),
                reference_reflections.flags.integrated_sum,
            )
        if (
            reference_crystal.get_space_group().type().number()
            != experiments.crystals()[0].get_space_group().type().number()
        ):
            raise Sorry("Space group of input does not match reference")
        try:
            reference_miller_set = filtered_arrays_from_experiments_reflections(
                reference_experiments, [reference_reflections]
            )[0]
        except ValueError:
            raise Sorry("No reflections remain after filtering the reference dataset")

        try:
            change_of_basis_op = change_of_basis_op_against_reference(
                experiments, reflections, reference_miller_set
            )
        except ValueError:
            raise Sorry("No reflections remain after filtering the test dataset")

    elif params.reference.file:
        wavelength = np.mean([expt.beam.get_wavelength() for expt in experiments])

        reference_miller_set = intensities_from_reference_file(
            params.reference.file, wavelength=wavelength
        )
        change_of_basis_op = change_of_basis_op_against_reference(
            experiments, reflections, reference_miller_set
        )

    elif len(experiments) and params.change_of_basis_op is Auto:
        if reference_crystal is not None:
            if len(experiments.crystals()) > 1:
                raise Sorry("Only one crystal can be processed at a time")
            from dials.algorithms.indexing.compare_orientation_matrices import (
                difference_rotation_matrix_axis_angle,
            )

            cryst = experiments.crystals()[0]
            R, axis, angle, change_of_basis_op = difference_rotation_matrix_axis_angle(
                cryst, reference_crystal
            )
            Rfmt = R.mathematica_form(format="%.3f", one_row_per_line=True)
            logger.info(
                "\n".join(
                    [
                        f"Change of basis op: {change_of_basis_op}",
                        "Rotation matrix to transform input crystal to reference::",
                        f"{Rfmt}",
                        f"Rotation of {angle:.3f} degrees",
                        f"about axis ({axis[0]:.3f}, {axis[1]:.3f}, {axis[2]:.3f})",
                    ]
                )
            )

        elif len(reflections):
            assert len(reflections) == 1

            # always re-map reflections to reciprocal space
            refl = reflections[0].deep_copy()
            refl.centroid_px_to_mm(experiments)
            refl.map_centroids_to_reciprocal_space(experiments)

            # index the reflection list using the input experiments list
            refl["id"] = flex.int(len(refl), -1)
            index = AssignIndicesGlobal(tolerance=0.2)
            index(refl, experiments)
            hkl_expt = refl["miller_index"]
            hkl_input = reflections[0]["miller_index"]

            change_of_basis_op = derive_change_of_basis_op(hkl_input, hkl_expt)

            # reset experiments list since we don't want to reindex this
            experiments = []

    else:
        change_of_basis_op = sgtbx.change_of_basis_op(params.change_of_basis_op)

    if len(experiments):
        space_group = params.space_group
        if space_group is not None:
            space_group = space_group.group()
        try:
            experiments = reindex_experiments(
                experiments, change_of_basis_op, space_group=space_group
            )
        except RuntimeError as e:
            # Only catch specific errors here
            if "Unsuitable value for rational rotation matrix." in str(e):
                original_message = str(e).split(":")[-1].strip()
                sys.exit(f"Error: {original_message} Is your change_of_basis_op valid?")
            elif "DXTBX_ASSERT(detail::is_r3_rotation_matrix(U)) failure" in str(e):
                sys.exit(
                    "Error: U is not a rotation matrix. Is your change_of_basis_op valid?"
                )
            raise

        logger.info(
            f"Saving reindexed experimental models to {params.output.experiments}"
        )
        experiments.as_file(params.output.experiments)

    if len(reflections):
        reflections = reindex_reflections(
            reflections, change_of_basis_op, params.hkl_offset
        )

        logger.info(f"Saving reindexed reflections to {params.output.reflections}")
        reflections.as_file(params.output.reflections)


if __name__ == "__main__":
    run()

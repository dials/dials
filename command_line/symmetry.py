from __future__ import division, absolute_import, print_function

import logging

logger = logging.getLogger("dials.command_line.symmetry")

import copy

from cctbx import sgtbx
import iotbx.phil

from dials.array_family import flex
from dials.util import log, Sorry
from dials.util.options import OptionParser, flatten_experiments, flatten_reflections
from dials.util.multi_dataset_handling import (
    assign_unique_identifiers,
    parse_multiple_datasets,
)
from dials.util.filter_reflections import filtered_arrays_from_experiments_reflections
from dials.algorithms.symmetry.determine_space_group import determine_space_group


phil_scope = iotbx.phil.parse(
    """\
d_min = Auto
  .type = float(value_min=0)

min_i_mean_over_sigma_mean = 4
  .type = float(value_min=0)

min_cc_half = 0.6
  .type = float(value_min=0, value_max=1)

batch = None
  .type = ints(value_min=0, size=2)

normalisation = kernel quasi ml_iso *ml_aniso
  .type = choice

lattice_group = None
  .type = space_group

verbosity = 0
  .type = int(value_min=0)
  .help = "The verbosity level"

seed = 230
  .type = int(value_min=0)

relative_length_tolerance = 0.05
  .type = float(value_min=0)

absolute_angle_tolerance = 2
  .type = float(value_min=0)

partiality_threshold = 0.99
  .type = float
  .help = "Use only reflections with a partiality above this threshold."

output {
  log = dials.symmetry.log
    .type = str
  debug_log = dials.symmetry.debug.log
    .type = str
  experiments = "symmetrized.expt"
    .type = path
  reflections = "symmetrized.refl"
    .type = path
  json = dials.symmetry.json
    .type = path
}

""",
    process_includes=True,
)


class symmetry(object):
    def __init__(self, experiments, reflections, params=None):
        if params is None:
            params = phil_scope.extract()
        self._params = params

        # transform models into miller arrays
        n_datasets = len(experiments)
        datasets = filtered_arrays_from_experiments_reflections(
            experiments,
            reflections,
            outlier_rejection_after_filter=True,
            partiality_threshold=params.partiality_threshold,
        )
        if len(datasets) != n_datasets:
            raise ValueError(
                """Some datasets have no reflection after prefiltering, please check
input data and filtering settings e.g partiality_threshold"""
            )

        result = determine_space_group(
            datasets,
            normalisation=self._params.normalisation,
            d_min=self._params.d_min,
            min_i_mean_over_sigma_mean=self._params.min_i_mean_over_sigma_mean,
            relative_length_tolerance=self._params.relative_length_tolerance,
            absolute_angle_tolerance=self._params.absolute_angle_tolerance,
        )
        logger.info(result)

        if params.output.json is not None:
            result.as_json(filename=params.output.json)

        self._export_experiments_reflections(experiments, reflections, result)

    def _export_experiments_reflections(self, experiments, reflections, result):
        from dxtbx.serialize import dump
        from rstbx.symmetry.constraints import parameter_reduction

        reindexed_experiments = copy.deepcopy(experiments)
        reindexed_reflections = flex.reflection_table()
        cb_op_inp_best = (
            result.best_solution.subgroup["cb_op_inp_best"] * result.cb_op_inp_min
        )
        best_subsym = result.best_solution.subgroup["best_subsym"]
        for i, expt in enumerate(reindexed_experiments):
            expt.crystal = expt.crystal.change_basis(result.cb_op_inp_min)
            expt.crystal.set_space_group(sgtbx.space_group("P 1"))
            expt.crystal = expt.crystal.change_basis(
                result.best_solution.subgroup["cb_op_inp_best"]
            )
            expt.crystal.set_space_group(
                best_subsym.space_group().build_derived_acentric_group()
            )
            S = parameter_reduction.symmetrize_reduce_enlarge(
                expt.crystal.get_space_group()
            )
            S.set_orientation(expt.crystal.get_B())
            S.symmetrize()
            expt.crystal.set_B(S.orientation.reciprocal_matrix())
            reindexed_refl = copy.deepcopy(reflections[i])
            reindexed_refl["miller_index"] = cb_op_inp_best.apply(
                reindexed_refl["miller_index"]
            )
            reindexed_reflections.extend(reindexed_refl)
        logger.info(
            "Saving reindexed experiments to %s" % self._params.output.experiments
        )
        dump.experiment_list(reindexed_experiments, self._params.output.experiments)
        logger.info(
            "Saving %s reindexed reflections to %s"
            % (len(reindexed_reflections), self._params.output.reflections)
        )
        reindexed_reflections.as_pickle(self._params.output.reflections)


help_message = """
This program implements the methods of
`POINTLESS <http://www.ccp4.ac.uk/html/pointless.html>`_ (
`Evans, P. (2006). Acta Cryst. D62, 72-82. <https://doi.org/10.1107/S0907444905036693>`_ and
`Evans, P. R. (2011). Acta Cryst. D67, 282-292. <https://doi.org/10.1107/S090744491003982X>`_)
for scoring and determination of Laue group symmetry.

The program takes as input a set of one or more integrated experiments and
reflections.

Examples::

  dials.symmetry models.expt observations.refl

"""


def run(args):
    usage = "dials.symmetry [options] models.expt observations.refl"

    parser = OptionParser(
        usage=usage,
        phil=phil_scope,
        read_reflections=True,
        read_experiments=True,
        check_format=False,
        epilog=help_message,
    )

    params, _, args = parser.parse_args(
        args=args, show_diff_phil=False, return_unhandled=True
    )

    # Configure the logging
    log.config(params.verbosity, info=params.output.log, debug=params.output.debug_log)

    from dials.util.version import dials_version

    logger.info(dials_version())

    # Log the diff phil
    diff_phil = parser.diff_phil.as_str()
    if diff_phil != "":
        logger.info("The following parameters have been modified:\n")
        logger.info(diff_phil)

    if params.seed is not None:
        import random

        flex.set_random_seed(params.seed)
        random.seed(params.seed)

    if not params.input.experiments or not params.input.reflections:
        parser.print_help()
        sys.exit()

    experiments = flatten_experiments(params.input.experiments)
    reflections = flatten_reflections(params.input.reflections)

    reflections = parse_multiple_datasets(reflections)
    if len(experiments) != len(reflections):
        raise Sorry(
            "Mismatched number of experiments and reflection tables found: %s & %s."
            % (len(experiments), len(reflections))
        )
    try:
        experiments, reflections = assign_unique_identifiers(experiments, reflections)
        symmetry(experiments, reflections, params=params)
    except ValueError as e:
        raise Sorry(e)


if __name__ == "__main__":
    import sys

    run(sys.argv[1:])

from __future__ import absolute_import, division, print_function

import iotbx.phil
from cctbx.array_family import flex

help_message = """

Computes the change of basis operator that minimises the difference between
two orientation matrices, and calculates the rotation matrix and Euler angles
that relate the two resulting orientation matrices (after transformation by
the calculated change of basis operator). Optionally calculates the angle
between given miller indices for the respective orientation matrices.

Examples::

  dials.compare_orientation_matrices models.expt

  dials.compare_orientation_matrices models_1.expt models_2.expt

  dials.compare_orientation_matrices models_1.expt models_2.expt hkl=1,0,0

"""


phil_scope = iotbx.phil.parse(
    """
hkl = None
  .type = ints(size=3)
  .multiple=True
comparison = *pairwise sequential
  .type = choice
space_group = None
  .type = space_group
"""
)


def run(args):
    from dials.util.options import OptionParser
    from dials.util.options import flatten_experiments
    import libtbx.load_env

    usage = "%s [options] models.expt" % libtbx.env.dispatcher_name

    parser = OptionParser(
        usage=usage,
        phil=phil_scope,
        read_experiments=True,
        check_format=False,
        epilog=help_message,
    )

    params, options = parser.parse_args(show_diff_phil=True)
    experiments = flatten_experiments(params.input.experiments)
    if len(experiments) <= 1:
        parser.print_help()
        return

    hkl = flex.miller_index(params.hkl)

    import dials.algorithms.indexing.compare_orientation_matrices

    crystals = []
    for experiment in experiments:
        crystal = experiment.crystal
        if params.space_group is not None:
            crystal.set_space_group(params.space_group.group())
        crystals.append(crystal)

    rmd = dials.algorithms.indexing.compare_orientation_matrices.rotation_matrix_differences(
        crystals, miller_indices=hkl, comparison=params.comparison
    )

    print(rmd)


if __name__ == "__main__":
    import sys

    run(sys.argv[1:])

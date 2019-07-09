# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1

from __future__ import absolute_import, division, print_function

import sys

import libtbx.load_env
import iotbx.phil
from dials.util import Sorry
from dials.util.options import OptionParser
from dials.util.options import flatten_experiments

help_message = """
"""

phil_scope = iotbx.phil.parse(
    """
threshold = 5000
  .type = float(value_min=0)
  .help = 'Threshold value for the clustering'
plot {
  show = False
    .type = bool
  name = 'cluster_unit_cell.png'
    .type = path
  log = False
    .type = bool
    .help = 'Display the dendrogram with a log scale'
}
"""
)


def run(args):
    usage = "%s [options] models.expt" % (libtbx.env.dispatcher_name)

    parser = OptionParser(
        usage=usage,
        phil=phil_scope,
        read_experiments=True,
        check_format=False,
        epilog=help_message,
    )

    params, options, args = parser.parse_args(
        show_diff_phil=True, return_unhandled=True
    )
    experiments = flatten_experiments(params.input.experiments)
    crystal_symmetries = []

    if len(experiments) == 0:
        if len(args):
            import os
            from iotbx.reflection_file_reader import any_reflection_file

            for arg in args:
                assert os.path.isfile(arg), arg
                reader = any_reflection_file(arg)
                assert reader.file_type() == "ccp4_mtz"
                arrays = reader.as_miller_arrays(merge_equivalents=False)
                crystal_symmetries.append(arrays[0].crystal_symmetry())
        else:
            parser.print_help()
            exit(0)
    else:
        from cctbx import crystal

        crystal_symmetries = [
            crystal.symmetry(
                unit_cell=expt.crystal.get_unit_cell(),
                space_group=expt.crystal.get_space_group(),
            )
            for expt in experiments
        ]

    do_cluster_analysis(crystal_symmetries, params)


def do_cluster_analysis(crystal_symmetries, params):

    try:
        from xfel.clustering.cluster import Cluster
        from xfel.clustering.cluster_groups import unit_cell_info
    except ImportError:
        raise Sorry("cluster_unit_cell requires xfel module but is not available")

    ucs = Cluster.from_crystal_symmetries(crystal_symmetries)

    if params.plot.show or params.plot.name is not None:
        if not params.plot.show:
            import matplotlib

            # http://matplotlib.org/faq/howto_faq.html#generate-images-without-having-a-window-appear
            matplotlib.use("Agg")  # use a non-interactive backend
        import matplotlib.pyplot as plt

        plt.figure("Andrews-Bernstein distance dendogram", figsize=(12, 8))
        ax = plt.gca()
        clusters, cluster_axes = ucs.ab_cluster(
            params.threshold,
            log=params.plot.log,
            ax=ax,
            write_file_lists=False,
            # schnell=_args.schnell,
            doplot=True,
        )
        print(unit_cell_info(clusters))
        plt.tight_layout()
        if params.plot.name is not None:
            plt.savefig(params.plot.name)
        if params.plot.show:
            plt.show()

    else:
        clusters, cluster_axes = ucs.ab_cluster(
            params.threshold,
            log=params.plot.log,
            write_file_lists=False,
            # schnell=_args.schnell,
            doplot=False,
        )
        print(unit_cell_info(clusters))

    return clusters


if __name__ == "__main__":
    run(sys.argv[1:])

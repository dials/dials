# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1


from __future__ import annotations

import json
import logging
import os

import iotbx.mtz
import iotbx.phil
from cctbx import crystal
from dxtbx.model import ExperimentList

import dials.util
from dials.algorithms.clustering.plots import scipy_dendrogram_to_plotly_json
from dials.algorithms.clustering.unit_cell import cluster_unit_cells
from dials.array_family import flex
from dials.util import log
from dials.util.multi_dataset_handling import (
    assign_unique_identifiers,
    parse_multiple_datasets,
)
from dials.util.options import ArgumentParser, reflections_and_experiments_from_files
from dials.util.version import dials_version

logger = logging.getLogger("dials.command_line.cluster_unit_cell")

help_message = """
"""

phil_scope = iotbx.phil.parse(
    """
threshold = 5000
  .type = float(value_min=0)
  .help = 'Threshold value for the clustering'
linkage = *single ward
  .type = choice
  .help = "The type of linkage to use for hierarchical clustering"
plot {
  show = False
    .type = bool
  name = 'cluster_unit_cell.png'
    .type = path
  log = False
    .type = bool
    .help = 'Display the dendrogram with a log scale'
}
output {
  clusters = False
    .type = bool
    .help = "If True, clusters will be split at the threshold value and a pair"
            "of output files will be created for each cluster"
  log = dials.cluster_unit_cell.log
    .type = str
  json = None
    .type = path
    .help = "Filename to output JSON summary of clustering results"
}
"""
)


@dials.util.show_mail_handle_errors()
def run(args=None):
    usage = "dials.cluster_unit_cell [options] models.expt"

    parser = ArgumentParser(
        usage=usage,
        phil=phil_scope,
        read_experiments=True,
        read_reflections=True,
        check_format=False,
        epilog=help_message,
    )

    params, _, args = parser.parse_args(
        args, show_diff_phil=False, return_unhandled=True
    )

    # Configure the logging
    log.config(logfile=params.output.log)
    logger.info(dials_version())

    # Log the diff phil
    diff_phil = parser.diff_phil.as_str()
    if diff_phil != "":
        logger.info("The following parameters have been modified:\n")
        logger.info(diff_phil)

    reflections, experiments = reflections_and_experiments_from_files(
        params.input.reflections, params.input.experiments
    )
    crystal_symmetries = []

    if len(experiments) == 0:
        if not args:
            parser.print_help()
            exit(0)
        for arg in args:
            assert os.path.isfile(arg), arg
            mtz_object = iotbx.mtz.object(file_name=arg)
            arrays = mtz_object.as_miller_arrays(
                merge_equivalents=False, anomalous=False
            )
            crystal_symmetries.append(arrays[0].crystal_symmetry())
    else:
        crystal_symmetries = [
            crystal.symmetry(
                unit_cell=expt.crystal.get_unit_cell(),
                space_group=expt.crystal.get_space_group(),
            )
            for expt in experiments
        ]

    if len(crystal_symmetries) <= 1:
        logger.info(f"Cannot cluster only {len(crystal_symmetries)} crystals, exiting.")
        exit(0)
    clusters = do_cluster_analysis(crystal_symmetries, params)

    if params.output.clusters:
        if len(experiments) == 0:
            logger.info("Clustering output can only be generated for input .expt files")
            return
        # Possibilities: either same number of experiments and reflection files,
        # or just one reflection file containing multiple sequences, or no
        # reflections given
        # Want a combined table to work on with experiment identifiers set.

        def _assign_and_return_joint(experiments, reflections):
            experiments, reflections = assign_unique_identifiers(
                experiments, reflections
            )
            joint_table = flex.reflection_table()
            for refls in reflections:
                joint_table.extend(refls)
            return joint_table

        if len(reflections) == 1:
            reflections = reflections[0]
            if len(set(reflections["id"])) != len(experiments):
                raise ValueError(
                    f"Mismatched number of reflection tables (f{len(set(reflections['id']))}) and experiments (f{len(experiments)})"
                )
            if not dict(reflections.experiment_identifiers()):
                reflections = reflections.split_by_experiment_id()
                reflections = _assign_and_return_joint(experiments, reflections)

        elif len(reflections) > 1:
            if not len(reflections) == len(experiments):
                reflections = parse_multiple_datasets(reflections)
                if len(reflections) != len(experiments):
                    raise ValueError(
                        f"Mismatched number of reflection tables (f{len(reflections)}) and experiments (f{len(experiments)})"
                    )
            reflections = _assign_and_return_joint(experiments, reflections)
        # else: no reflections given, continue and just split experiments

        clusters.sort(key=len, reverse=True)
        for cluster in clusters:
            sub_expt = ExperimentList([experiments[i] for i in cluster.lattice_ids])
            expt_filename = cluster.name + ".expt"
            logger.info(
                f"Saving {len(sub_expt)} lattices from {cluster.name} to {expt_filename}"
            )
            sub_expt.as_file(expt_filename)
            if reflections:
                identifiers = sub_expt.identifiers()
                sub_refl = reflections.select_on_experiment_identifiers(identifiers)
                # renumber the ids to go from 0->n-1
                sub_refl.reset_ids()
                refl_filename = cluster.name + ".refl"
                logger.info(
                    f"Saving reflections from {cluster.name} to {refl_filename}"
                )
                sub_refl.as_file(refl_filename)


def do_cluster_analysis(crystal_symmetries, params):
    lattice_ids = list(range(len(crystal_symmetries)))

    if params.plot.show or params.plot.name:
        if not params.plot.show:
            import matplotlib

            # http://matplotlib.org/faq/howto_faq.html#generate-images-without-having-a-window-appear
            matplotlib.use("Agg")  # use a non-interactive backend
        import matplotlib.pyplot as plt

        plt.figure("Andrews-Bernstein distance dendogram", figsize=(12, 8))
        ax = plt.gca()
        no_plot = False
    else:
        ax = None
        no_plot = True

    clustering = cluster_unit_cells(
        crystal_symmetries,
        lattice_ids=lattice_ids,
        threshold=params.threshold,
        ax=ax,
        no_plot=no_plot,
        linkage=params.linkage,
    )
    logger.info(clustering)

    if params.plot.show or params.plot.name:
        if params.plot.log:
            ax.set_yscale("symlog", linthresh=1)
        else:
            ax.set_ylim(-ax.get_ylim()[1] / 100, ax.get_ylim()[1])

        plt.tight_layout()
        if params.plot.name:
            plt.savefig(params.plot.name)
        if params.plot.show:
            plt.show()

    if params.output.json:
        dendrogram_json = scipy_dendrogram_to_plotly_json(
            clustering.dendrogram,
            title="Unit cell clustering",
            xtitle="Dataset",
            ytitle="Distance (Å<sup>2</sup>)",
            help="""\
    The results of single-linkage hierarchical clustering on the unit cell parameters using
    the Andrews–Bernstein NCDist distance metric (Andrews & Bernstein, 2014). The height at
    which two clusters are merged in the dendrogram is a measure of the similarity between
    the unit cells in each cluster. A larger separation between two clusters may be
    indicative of a higher degree of non-isomorphism between the clusters. Conversely, a
    small separation between two clusters suggests that their unit cell parameters are
    relatively isomorphous.
    """,
        )

        with open(params.output.json, "w") as f:
            json.dump(dendrogram_json, f, indent=2)
        logger.info(f"Wrote clustering dendrogram JSON to {params.output.json}")

    return clustering.clusters


if __name__ == "__main__":
    run()

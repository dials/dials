from __future__ import annotations

import logging
import os
import sys
from collections.abc import Sequence
from typing import Iterator, List, Optional, TypeVar

from dxtbx.model.experiment_list import ExperimentList
from libtbx import phil
from libtbx.phil import parse

import dials.util
from dials.array_family import flex
from dials.util import log
from dials.util.combine_experiments import CombineWithReference  # noqa
from dials.util.combine_experiments import (
    combine_experiments,
    combine_experiments_no_reflections,
    do_unit_cell_clustering,
)
from dials.util.options import ArgumentParser
from dials.util.version import dials_version

T = TypeVar("T")
logger = logging.getLogger(__name__)

help_message = """

Utility script to combine multiple reflections and experiments files into
one multi-experiment reflections and one experiments file. Experiments are
matched to reflections in the order they are provided as input.

Reference models can be chosen from any of the input experiments files. These
will replace all other models of that type in the output experiments file.
This is useful, for example, for combining multiple experiments that should
differ only in their crystal models. No checks are made to ensure that a
reference model is a good replacement model.

Although only one reference model of each type is allowed, more complex
combinations of experiments can be created by repeat runs.

Examples::

  dials.combine_experiments experiments_0.expt experiments_1.expt \\
    reflections_0.refl reflections_1.refl \\
    reference_from_experiment.beam=0 \\
    reference_from_experiment.detector=0
"""

# The phil scope
phil_scope = parse(
    """
  output {
    log = dials.combine_experiments.log
      .type = str
  }
  reference_from_experiment{
    beam = None
      .help = "Take beam model from this experiment to overwrite all other"
              "beam models in the combined experiments"
      .type = int(value_min=0)

    scan = None
      .help = "Take scan model from this experiment to overwrite all other"
              "scan models in the combined experiments"
      .type = int(value_min=0)

    crystal = None
      .help = "Take crystal model from this experiment to overwrite all"
              "other crystal models in the combined experiments"
      .type = int(value_min=0)

    goniometer = None
      .help = "Take goniometer model from this experiment to overwrite all"
              "other goniometer models in the combined experiments"
      .type = int(value_min=0)

    detector = None
      .help = "Take detector model from this experiment to overwrite all"
              "other detector models in the combined experiments"
      .type = int(value_min=0)

    average_detector = False
      .help = "Create an average detector model from all the input detector"
              "models and use it as the reference. Not compatible with"
              "reference_from_experiment.detector"
      .type = bool

    compare_models = True
      .help = "Whether to compare a model with the reference model before"
              "replacing it. If the comparison falls outside the tolerance,"
              "the combination will not be allowed. Disable comparison to force"
              "overwriting of models with the reference"
      .type = bool

    average_hierarchy_level = None
      .help = "For hierarchical detectors, optionally provide a single level"
              "to do averaging at."
      .type = int(value_min=0)

    include scope dials.util.options.tolerance_phil_scope

  }

  clustering {
    use = False
      .type = bool
      .help = "Separate experiments into subsets using the clustering"
              "toolkit. One json per cluster will be saved."

    dendrogram = False
      .type = bool
      .help = "Display dendrogram of the clustering results. Should not"
              "be used with parallel processing."

    threshold = 1000
      .type = int
      .help = "Threshold used in the dendrogram to separate into clusters."

    max_clusters = None
      .type = int
      .help = "Maximum number of clusters to save as jsons."

    exclude_single_crystal_clusters = True
      .type = bool
      .help = "Don't produce a 'cluster' containing only one crystal."

  }

  output {
    experiments_filename = combined.expt
      .type = str
      .help = "The filename for combined experimental models"

    reflections_filename = combined.refl
      .type = str
      .help = "The filename for combined reflections"

    n_subset = None
      .type = int
      .help = "If not None, keep a subset of size n_subset when"
              "saving the combined experiments"

    n_subset_method = *random n_refl significance_filter
      .type = choice
      .help = "Algorithm to be used for choosing the n_subset images/"
              "experiments for refinement.  n_refl chooses the set with the"
              "largest numbers of reflections listed in the pickle files"
              "significance filter used to select n_subset images based on"
              "I/sig(I) cutoff"

    n_refl_panel_list = None
      .type = ints
      .help = "If n_subset_method is n_refl, specify which panels to search"
              "on."

    max_batch_size = None
      .type = int
      .expert_level = 2
      .help = "If not None, split the resultant combined set of experiments"
              "into separate files, each at most max_batch_size number of"
              "experiments. Example, if there were 5500 experiments and"
              "max_batch_size is 1000, 6 experiment lists will be created,"
              "of sizes 917, 917, 917, 917, 916, 916"

    delete_shoeboxes = False
      .type = bool
      .expert_level = 2
      .help = "If true, delete shoeboxes from reflection tables while comb-"
              "ining them to save on memory."

    min_reflections_per_experiment = None
      .type = int
      .expert_level = 2
      .help = "If not None, throw out any experiment with fewer than this"
              "many reflections"

    max_reflections_per_experiment = None
      .type = int
      .expert_level = 2
      .help = "If not None, throw out any experiment with more than this"
              "many reflections"

    include scope dials.algorithms.integration.stills_significance_filter.phil_scope
  }
""",
    process_includes=True,
)


def _split_equal_parts_of_length(a: Sequence[T], n: int) -> Iterator[Sequence[T]]:
    k, m = divmod(len(a), n)
    return (a[i * k + min(i, m) : (i + 1) * k + min(i + 1, m)] for i in range(n))


def save_combined_experiments(
    experiments,
    clustering: phil.scope_extract,
    reflections: Optional[flex.reflection_table] = None,
    max_batch_size: int = None,
    experiments_filename="combined.expt",
    reflections_filename="combined.refl",
):
    output_experiments_list: List[ExperimentList] = []
    output_reflections_list: List[flex.reflection_table] = []
    expt_names_list: List[str] = []
    refl_names_list: List[str] = []

    # cluster the resulting experiments if requested
    if clustering.use and len(experiments) > 1:
        if clustering.max_clusters == 0:
            sys.exit("Error: max_clusters must be None or >0")
        clustered = do_unit_cell_clustering(
            experiments,
            dendrogram=clustering.dendrogram,
            threshold=clustering.threshold,
        )
        n_clusters = len(clustered)
        clusters = sorted(clustered.clusters, key=len, reverse=True)
        if clustering.max_clusters is not None:
            clusters = clusters[: clustering.max_clusters]
        if clustering.exclude_single_crystal_clusters:
            clusters = [c for c in clusters if len(c) > 1]
        output_experiments_list = [
            ExperimentList(
                [expt for expt in experiments if expt.identifier in cluster.lattice_ids]
            )
            for cluster in clusters
        ]
        if reflections:
            output_reflections_list = [
                reflections.select_on_experiment_identifiers(cluster.lattice_ids)
                for cluster in clusters
            ]
            for table in output_reflections_list:
                table.reset_ids()
        for i_cluster, _ in enumerate(output_experiments_list):
            exp_name = os.path.splitext(experiments_filename)[0] + (
                "_cluster%d.expt" % (n_clusters - i_cluster)
            )
            expt_names_list.append(exp_name)
            if reflections:
                refl_name = os.path.splitext(reflections_filename)[0] + (
                    "_cluster%d.refl" % (n_clusters - i_cluster)
                )
                refl_names_list.append(refl_name)

    else:
        output_experiments_list = [experiments]
        expt_names_list = [experiments_filename]
        if reflections:
            output_reflections_list = [reflections]
            refl_names_list = [reflections_filename]

    if reflections:
        _save_experiments_and_reflections(
            output_experiments_list,
            expt_names_list,
            output_reflections_list,
            refl_names_list,
            max_batch_size,
        )
    else:
        _save_only_experiments(output_experiments_list, expt_names_list, max_batch_size)


def _save_only_experiments(
    output_experiments_list, expt_names_list, max_batch_size=None
):
    for elist, ename in zip(output_experiments_list, expt_names_list):
        if max_batch_size is None:
            logger.info(f"Saving combined experiments to {ename}")
            elist.as_file(ename)
        else:
            for i, indices in enumerate(
                _split_equal_parts_of_length(
                    list(range(len(elist))), (len(elist) // max_batch_size) + 1
                )
            ):
                batch_expts = ExperimentList()
                for sub_id, sub_idx in enumerate(indices):
                    batch_expts.append(elist[sub_idx])
                exp_filename = os.path.splitext(ename)[0] + "_%03d.expt" % i
                logger.info(f"Saving combined experiments to {exp_filename}")
                batch_expts.as_file(exp_filename)


def _save_experiments_and_reflections(
    output_experiments_list,
    expt_names_list,
    output_reflections_list,
    refl_names_list,
    max_batch_size=None,
):
    for elist, table, ename, rname in zip(
        output_experiments_list,
        output_reflections_list,
        expt_names_list,
        refl_names_list,
    ):
        if max_batch_size is None:
            logger.info(f"Saving combined experiments to {ename}")
            elist.as_file(ename)
            logger.info(f"Saving combined reflections to {rname}")
            table.as_file(rname)
        else:
            for i, indices in enumerate(
                _split_equal_parts_of_length(
                    list(range(len(elist))), (len(elist) // max_batch_size) + 1
                )
            ):
                batch_expts = ExperimentList()
                batch_refls = flex.reflection_table()
                if table.experiment_identifiers().keys():
                    for sub_idx in indices:
                        batch_expts.append(elist[sub_idx])
                    batch_refls = table.select(batch_expts)
                    batch_refls.reset_ids()
                else:
                    for sub_id, sub_idx in enumerate(indices):
                        batch_expts.append(elist[sub_idx])
                        sub_refls = table.select(table["id"] == sub_idx)
                        sub_refls["id"] = flex.int(len(sub_refls), sub_id)
                        batch_refls.extend(sub_refls)
                exp_filename = os.path.splitext(ename)[0] + "_%03d.expt" % i
                ref_filename = os.path.splitext(rname)[0] + "_%03d.refl" % i
                logger.info(f"Saving combined experiments to {exp_filename}")
                batch_expts.as_file(exp_filename)
                logger.info(f"Saving combined reflections to {ref_filename}")
                batch_refls.as_file(ref_filename)


@dials.util.show_mail_handle_errors()
def run(args=None) -> None:
    usage = (
        "usage: dials.combine_experiments [options] [param.phil] "
        "experiments1.expt experiments2.expt reflections1.refl "
        "reflections2.refl..."
    )

    # Create the parser
    parser = ArgumentParser(
        usage=usage,
        phil=phil_scope,
        read_reflections=True,
        read_experiments=True,
        check_format=False,
        epilog=help_message,
    )
    params, options = parser.parse_args(args, show_diff_phil=True)

    # Configure the logging
    log.config(verbosity=options.verbose, logfile=params.output.log)
    logger.info(dials_version())

    if not params.input.experiments:
        logger.info("No Experiments found in the input")
        parser.print_help()
        return
    if not params.input.reflections:
        if (
            params.output.n_subset is not None
            and params.output.n_subset_method != "random"
        ):
            logger.info(
                """No reflection data found in the input.
Reflection tables are needed if n_subset_method != random and n_subset is not None"""
            )
            parser.print_help()
            return
    else:
        if len(params.input.reflections) != len(params.input.experiments):
            sys.exit(
                "The number of input reflections files does not match the "
                "number of input experiments"
            )

    # we want a list of experiment lists and list of experiment tables (one per input file)
    experiment_lists: List[ExperimentList] = [
        ExperimentList(o.data) for o in params.input.experiments
    ]
    reflection_tables: List[flex.reflection_table] = [
        o.data for o in params.input.reflections
    ]
    if reflection_tables:
        expts, refls = combine_experiments(params, experiment_lists, reflection_tables)
    else:
        expts = combine_experiments_no_reflections(params, experiment_lists)
        refls = None
    save_combined_experiments(
        expts,
        params.clustering,
        refls,
        max_batch_size=params.output.max_batch_size,
        experiments_filename=params.output.experiments_filename,
        reflections_filename=params.output.reflections_filename,
    )


if __name__ == "__main__":
    run()

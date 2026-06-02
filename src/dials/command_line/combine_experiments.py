from __future__ import annotations

import logging
import os
import sys
from collections.abc import Iterator, Sequence
from typing import TypeVar

from dxtbx.model.experiment_list import ExperimentList
from libtbx import phil
from libtbx.phil import parse

import dials.util
from dials.array_family import flex
from dials.util import log
from dials.util.combine_experiments import (
    CombineWithReference,  # noqa
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

    sort_by_imageset_path_and_image_index = False
      .type = bool
      .expert_level = 2
      .help = "If True, sort the experiments and reflections first by path"
              "then by image number (for composite files like HDF5)"

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
    reflections: flex.reflection_table | None = None,
    max_batch_size: int = None,
    experiments_filename="combined.expt",
    reflections_filename="combined.refl",
):
    output_experiments_list: list[ExperimentList] = []
    output_reflections_list: list[flex.reflection_table] = []
    expt_names_list: list[str] = []
    refl_names_list: list[str] = []

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


def _consolidate_stills_imagesets(experiments, reflections=None):
    """Merge experiments from the same source file into one shared ImageSequence.

    When N per-worker .expt files are combined, each file contributes its own
    ImageSequence Python object even when all files reference the same source
    data file. This causes to_dict() to write N imageset entries. Group by
    source path and re-create one shared ImageSequence per unique path.

    Preserves experiment order (so reflection-table 'id' stays valid).
    Does not re-open source files.
    """
    from collections import defaultdict

    from dxtbx.imageset import ImageSequence
    from dxtbx.model import Scan
    from dxtbx.model.experiment_list import Experiment

    if not experiments:
        return experiments, reflections
    first_iset = experiments[0].imageset
    if not isinstance(first_iset, ImageSequence):
        return experiments, reflections
    if first_iset.get_scan() is None or not first_iset.get_scan().is_still():
        return experiments, reflections

    by_path = defaultdict(list)
    for expt in experiments:
        by_path[expt.imageset.paths()[0]].append(expt)

    def _frame(expt):
        # 0-based frame this experiment owns within its source file. Derive it
        # from the per-experiment scan (the same key the to_dict consolidation
        # uses to assign scan_point and that the reader zips against
        # single_file_indices), NOT from imageset.indices(): a shared
        # ImageSequence reports all its frames via indices(), but after
        # min/max_reflections_per_experiment / n_subset filtering only a subset
        # of experiments survives. image_range[0] - 1 is the value that
        # round-trips back to this experiment's scan on read.
        if expt.scan is not None:
            return expt.scan.get_image_range()[0] - 1
        return expt.imageset.indices()[0]

    def _needs_rebuild(grp):
        # Cross-worker merge: distinct imageset objects for one source path must
        # be collapsed to one shared ImageSequence.
        if len({id(e.imageset) for e in grp}) > 1:
            return True
        # A single-frame imageset (e.g. multi-file CBF stills, one .cbf per
        # shot) has nothing to prune: indices() == [0] and its scan image_range
        # is the absolute CBF frame number from the filename, not a 0-based
        # index into a shared file, so the indices()-vs-_frame() comparison
        # below would always (wrongly) mismatch. Only a genuine multi-frame
        # composite can over-advertise frames. (The cross-worker branch above
        # also never fires here: one experiment per file => one path group per
        # imageset, never a shared path with distinct imageset objects.)
        iset = grp[0].imageset
        if len(iset) <= 1:
            return False
        # Prune-after-filter: a single shared composite that still advertises
        # frames with no surviving experiment must be pruned, or the consolidated
        # scan's single_file_indices over-counts the scan_points on read.
        return set(iset.indices()) != {_frame(e) for e in grp}

    # Nothing to do if no path group needs merging or pruning (the common
    # unfiltered single-imageset combine).
    if not any(_needs_rebuild(grp) for grp in by_path.values()):
        return experiments, reflections

    old_imagesets = list(experiments.imagesets())

    from dxtbx.format.Registry import get_format_class_for_file

    path_to_iset = {}
    for path, grp in by_path.items():
        # Use only the surviving experiments' frames (sparse), matching the
        # sparse layout produced by stills_process._rebuild_shared_imageset_output
        # so downstream consumers (image_viewer's frame-list lookup, etc.) see a
        # consistent imageset structure regardless of which tool produced it.
        sorted_frames = sorted({_frame(e) for e in grp})
        if len({id(e.imageset) for e in grp}) > 1:
            # Cross-worker merge: each worker's imageset data() is sized for that
            # worker's frame subset, so re-open the source file to get a
            # full-file ImageSetData that can accommodate the union of frames.
            data = get_format_class_for_file(path).get_imageset([path]).data()
        else:
            # Prune: the single shared imageset's data() already covers all its
            # own frames and we are only removing some, so reuse it. Re-opening
            # is unnecessary (and would fail for synthetic readers whose path
            # does not exist on disk).
            data = grp[0].imageset.data()
        path_to_iset[path] = ImageSequence(
            data,
            flex.size_t(sorted_frames),
            grp[0].beam,
            grp[0].detector,
            None,
            Scan((1, len(sorted_frames)), (0.0, 0.0)),
        )

    new_experiments = ExperimentList()
    for expt in experiments:
        new_experiments.append(
            Experiment(
                imageset=path_to_iset[expt.imageset.paths()[0]],
                beam=expt.beam,
                detector=expt.detector,
                goniometer=expt.goniometer,
                scan=expt.scan,
                crystal=expt.crystal,
                profile=expt.profile,
                scaling_model=expt.scaling_model,
                identifier=expt.identifier,
            )
        )

    if reflections is not None and "imageset_id" in reflections:
        new_imagesets = list(new_experiments.imagesets())
        old_to_new = {
            old_imagesets.index(old_e.imageset): new_imagesets.index(new_e.imageset)
            for old_e, new_e in zip(experiments, new_experiments)
        }
        new_iset_id = reflections["imageset_id"].deep_copy()
        for old_id, new_id in old_to_new.items():
            new_iset_id.set_selected(reflections["imageset_id"] == old_id, new_id)
        reflections["imageset_id"] = new_iset_id

    return new_experiments, reflections


def _sort_experiments_and_reflections(expts, refls):
    """Sort experiments and reflections by (source path, frame index).

    Works for both per-frame imagesets (imageset.indices()[0]) and shared
    ImageSequences produced by _consolidate_stills_imagesets (uses the
    per-experiment scan's image_range).  Remaps reflection-table 'id',
    'experiment_identifiers', and 'imageset_id' to stay consistent with
    the new experiment order.
    """
    from dxtbx.imageset import ImageSequence

    def _sort_key(expt):
        iset = expt.imageset
        path = iset.paths()[0]
        if isinstance(iset, ImageSequence) and len(iset) > 1:
            fi = expt.scan.get_image_range()[0] if expt.scan is not None else 0
        else:
            fi = iset.indices()[0]
        return (path, fi)

    perm = sorted(range(len(expts)), key=lambda i: _sort_key(expts[i]))
    if perm == list(range(len(expts))):
        return expts, refls

    old_imagesets = list(expts.imagesets())
    sorted_expts = ExperimentList([expts[perm[i]] for i in range(len(perm))])
    new_imagesets = list(sorted_expts.imagesets())
    iset_remap = {
        old_id: new_id
        for new_id, iset in enumerate(new_imagesets)
        for old_id, old_iset in enumerate(old_imagesets)
        if iset is old_iset
    }

    if refls:
        old_idents = dict(refls.experiment_identifiers())
        new_blocks = []
        for new_id, old_id in enumerate(perm):
            block = refls.select(refls["id"] == old_id)
            if len(block):
                block["id"] = flex.int(len(block), new_id)
                new_blocks.append(block)
        if new_blocks:
            sorted_refls = flex.reflection_table.concat(new_blocks)
            idents = sorted_refls.experiment_identifiers()
            for k in list(idents.keys()):
                del idents[k]
            for new_id, old_id in enumerate(perm):
                if old_id in old_idents:
                    idents[new_id] = old_idents[old_id]
            if "imageset_id" in sorted_refls and any(
                o != n for o, n in iset_remap.items()
            ):
                old_iset_id = sorted_refls["imageset_id"].deep_copy()
                new_iset_id = old_iset_id.deep_copy()
                for old_id, new_id in iset_remap.items():
                    if old_id != new_id:
                        new_iset_id.set_selected(old_iset_id == old_id, new_id)
                sorted_refls["imageset_id"] = new_iset_id
        else:
            sorted_refls = refls
    else:
        sorted_refls = refls
    return sorted_expts, sorted_refls


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
    experiment_lists: list[ExperimentList] = [
        ExperimentList(o.data) for o in params.input.experiments
    ]
    reflection_tables: list[flex.reflection_table] = [
        o.data for o in params.input.reflections
    ]
    if reflection_tables:
        expts, refls = combine_experiments(params, experiment_lists, reflection_tables)
    else:
        expts = combine_experiments_no_reflections(params, experiment_lists)
        refls = None
    expts, refls = _consolidate_stills_imagesets(expts, refls)
    _is_stills = bool(expts) and expts[0].is_still()
    if _is_stills or params.output.sort_by_imageset_path_and_image_index:
        expts, refls = _sort_experiments_and_reflections(expts, refls)
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

from __future__ import annotations

import logging
import os
import random
import sys
from typing import List, Tuple

import dxtbx.model
from dxtbx.model.experiment_list import ExperimentList
from libtbx.phil import parse
from scitbx import matrix

import dials.util
from dials.algorithms.integration.stills_significance_filter import SignificanceFilter
from dials.array_family import flex
from dials.util import log, tabulate
from dials.util.combine_experiments import (
    CombineWithReference,
    ComparisonError,
    _split_equal_parts_of_length,
    do_unit_cell_clustering,
    find_experiment_in,
)
from dials.util.options import ArgumentParser, flatten_experiments
from dials.util.version import dials_version

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


def average_detectors(target, panelgroups, depth, average_hierarchy_level=None):
    # Recursive function to do the averaging

    if average_hierarchy_level is None or depth == average_hierarchy_level:
        n = len(panelgroups)
        sum_fast = matrix.col((0.0, 0.0, 0.0))
        sum_slow = matrix.col((0.0, 0.0, 0.0))
        sum_ori = matrix.col((0.0, 0.0, 0.0))

        # Average the d matrix vectors
        for pg in panelgroups:
            sum_fast += matrix.col(pg.get_local_fast_axis())
            sum_slow += matrix.col(pg.get_local_slow_axis())
            sum_ori += matrix.col(pg.get_local_origin())
        sum_fast /= n
        sum_slow /= n
        sum_ori /= n

        # Re-orthagonalize the slow and the fast vectors by rotating around the cross product
        c = sum_fast.cross(sum_slow)
        a = sum_fast.angle(sum_slow, deg=True) / 2
        sum_fast = sum_fast.rotate_around_origin(c, a - 45, deg=True)
        sum_slow = sum_slow.rotate_around_origin(c, -(a - 45), deg=True)

        target.set_local_frame(sum_fast, sum_slow, sum_ori)

    if target.is_group():
        # Recurse
        for i, target_pg in enumerate(target):
            average_detectors(
                target_pg,
                [pg[i] for pg in panelgroups],
                depth + 1,
                average_hierarchy_level,
            )


def parse_ref_models(
    flat_exps: ExperimentList, reference_from_experiment: phil_scope
) -> Tuple[
    None | dxtbx.model.beam,
    None | dxtbx.model.goniometer,
    None | dxtbx.model.scan,
    None | dxtbx.model.crystal,
    None | dxtbx.model.detector,
]:
    ref_beam = reference_from_experiment.beam
    ref_goniometer = reference_from_experiment.goniometer
    ref_scan = reference_from_experiment.scan
    ref_crystal = reference_from_experiment.crystal
    ref_detector = reference_from_experiment.detector

    if ref_beam is not None:
        try:
            ref_beam = flat_exps[ref_beam].beam
        except IndexError:
            sys.exit(f"{ref_beam} is not a valid experiment ID")

    if ref_goniometer is not None:
        try:
            ref_goniometer = flat_exps[ref_goniometer].goniometer
        except IndexError:
            sys.exit(f"{ref_goniometer} is not a valid experiment ID")

    if ref_scan is not None:
        try:
            ref_scan = flat_exps[ref_scan].scan
        except IndexError:
            sys.exit(f"{ref_scan} is not a valid experiment ID")

    if ref_crystal is not None:
        try:
            ref_crystal = flat_exps[ref_crystal].crystal
        except IndexError:
            sys.exit(f"{ref_crystal} is not a valid experiment ID")

    if ref_detector is not None:
        assert not reference_from_experiment.average_detector
        try:
            ref_detector = flat_exps[ref_detector].detector
        except IndexError:
            sys.exit(f"{ref_detector} is not a valid experiment ID")
    elif reference_from_experiment.average_detector:
        # Average all of the detectors together
        ref_detector = flat_exps[0].detector
        average_detectors(
            ref_detector.hierarchy(),
            [e.detector.hierarchy() for e in flat_exps],
            0,
            reference_from_experiment.average_hierarchy_level,
        )
    return (ref_beam, ref_goniometer, ref_scan, ref_crystal, ref_detector)


def run_with_preparsed(params, flat_exps):
    """Run combine_experiments, but allow passing in of parameters"""

    ref_beam, ref_goniometer, ref_scan, ref_crystal, ref_detector = parse_ref_models(
        flat_exps, params.reference_from_experiment
    )

    combine = CombineWithReference(
        beam=ref_beam,
        goniometer=ref_goniometer,
        scan=ref_scan,
        crystal=ref_crystal,
        detector=ref_detector,
        params=params,
    )

    # set up global experiments and reflections lists
    reflections = flex.reflection_table()
    global_id = 0
    skipped_expts_min_refl = 0
    skipped_expts_max_refl = 0
    experiments = ExperimentList()

    # loop through the input, building up the global lists
    nrefs_per_exp = []
    for ref_wrapper, exp_wrapper in zip(
        params.input.reflections, params.input.experiments
    ):
        refs = ref_wrapper.data
        exps = exp_wrapper.data

        # Record initial mapping of ids for updating later.
        ids_map = dict(refs.experiment_identifiers())

        # Keep track of mapping of imageset_ids old->new within this experimentlist
        imageset_result_map = {}

        for k in refs.experiment_identifiers().keys():
            del refs.experiment_identifiers()[k]
        for i, exp in enumerate(exps):
            sel = refs["id"] == i
            sub_ref = refs.select(sel)
            n_sub_ref = len(sub_ref)
            if (
                params.output.min_reflections_per_experiment is not None
                and n_sub_ref < params.output.min_reflections_per_experiment
            ):
                skipped_expts_min_refl += 1
                continue
            if (
                params.output.max_reflections_per_experiment is not None
                and n_sub_ref > params.output.max_reflections_per_experiment
            ):
                skipped_expts_max_refl += 1
                continue

            nrefs_per_exp.append(n_sub_ref)
            sub_ref["id"] = flex.int(len(sub_ref), global_id)

            # now update identifiers if set.
            if i in ids_map:
                sub_ref.experiment_identifiers()[global_id] = ids_map[i]
            if params.output.delete_shoeboxes and "shoebox" in sub_ref:
                del sub_ref["shoebox"]

            try:
                experiments.append(combine(exp))
            except ComparisonError as e:
                # When we failed tolerance checks, give a useful error message
                (path, index) = find_experiment_in(exp, params.input.experiments)
                sys.exit(
                    "Model didn't match reference within required tolerance for experiment {} in {}:"
                    "\n{}\nAdjust tolerances or set compare_models=False to ignore differences.".format(
                        index, path, str(e)
                    )
                )

            # Rewrite imageset_id, if the experiment has and imageset
            if exp.imageset and "imageset_id" in sub_ref:
                # Get the index of the imageset for this experiment and record how it changed
                new_imageset_id = experiments.imagesets().index(
                    experiments[-1].imageset
                )
                old_imageset_id = exps.imagesets().index(exp.imageset)
                imageset_result_map[old_imageset_id] = new_imageset_id

                # Check for invalid(?) imageset_id indices... and leave if they are wrong
                if len(set(sub_ref["imageset_id"])) != 1:
                    logger.warning(
                        "Warning: Experiment %d reflections appear to have come from multiple imagesets - output may be incorrect",
                        i,
                    )
                else:
                    sub_ref["imageset_id"] = flex.int(len(sub_ref), new_imageset_id)

            reflections.extend(sub_ref)

            global_id += 1

        # Include unindexed reflections, if we can safely remap their imagesets
        if "imageset_id" in reflections:
            unindexed_refs = refs.select(refs["id"] == -1)
            for old_id in set(unindexed_refs["imageset_id"]):
                subs = unindexed_refs.select(unindexed_refs["imageset_id"] == old_id)
                subs["imageset_id"] = flex.int(len(subs), imageset_result_map[old_id])
                reflections.extend(subs)

    if (
        params.output.min_reflections_per_experiment is not None
        and skipped_expts_min_refl > 0
    ):
        print(
            "Removed {} experiments with fewer than {} reflections".format(
                skipped_expts_min_refl, params.output.min_reflections_per_experiment
            )
        )
    if (
        params.output.max_reflections_per_experiment is not None
        and skipped_expts_max_refl > 0
    ):
        print(
            "Removed {} experiments with more than {} reflections".format(
                skipped_expts_max_refl, params.output.max_reflections_per_experiment
            )
        )

    # print number of reflections per experiment

    header = ["Experiment", "Number of reflections"]
    rows = [(str(i), str(n)) for (i, n) in enumerate(nrefs_per_exp)]
    print(tabulate(rows, header))

    # save a random subset if requested
    if params.output.n_subset is not None and len(experiments) > params.output.n_subset:
        subset_exp = ExperimentList()
        subset_refls = flex.reflection_table()
        if params.output.n_subset_method == "random":
            n_picked = 0
            indices = list(range(len(experiments)))
            if reflections.experiment_identifiers().keys():
                indices_to_sel = []
                while n_picked < params.output.n_subset:
                    idx = indices.pop(random.randint(0, len(indices) - 1))
                    indices_to_sel.append(idx)
                    n_picked += 1
                # make sure select in order.
                for idx in sorted(indices_to_sel):
                    subset_exp.append(experiments[idx])
                subset_refls = reflections.select(subset_exp)
                subset_refls.reset_ids()
            else:
                while n_picked < params.output.n_subset:
                    idx = indices.pop(random.randint(0, len(indices) - 1))
                    subset_exp.append(experiments[idx])
                    refls = reflections.select(reflections["id"] == idx)
                    refls["id"] = flex.int(len(refls), n_picked)
                    subset_refls.extend(refls)
                    n_picked += 1
            print(
                "Selecting a random subset of {} experiments out of {} total.".format(
                    params.output.n_subset, len(experiments)
                )
            )
        elif params.output.n_subset_method == "n_refl":
            if params.output.n_refl_panel_list is None:
                refls_subset = reflections
            else:
                sel = flex.bool(len(reflections), False)
                for p in params.output.n_refl_panel_list:
                    sel |= reflections["panel"] == p
                refls_subset = reflections.select(sel)
            refl_counts = flex.int()
            for expt_id in range(len(experiments)):
                refl_counts.append((refls_subset["id"] == expt_id).count(True))
            sort_order = flex.sort_permutation(refl_counts, reverse=True)
            if reflections.experiment_identifiers().keys():
                for idx in sorted(sort_order[: params.output.n_subset]):
                    subset_exp.append(experiments[idx])
                subset_refls = reflections.select(subset_exp)
                subset_refls.reset_ids()
            else:
                for expt_id, idx in enumerate(sort_order[: params.output.n_subset]):
                    subset_exp.append(experiments[idx])
                    refls = reflections.select(reflections["id"] == idx)
                    refls["id"] = flex.int(len(refls), expt_id)
                    subset_refls.extend(refls)
            print(
                "Selecting a subset of {} experiments with highest number of reflections out of {} total.".format(
                    params.output.n_subset, len(experiments)
                )
            )

        elif params.output.n_subset_method == "significance_filter":
            params.output.significance_filter.enable = True
            sig_filter = SignificanceFilter(params.output)
            refls_subset = sig_filter(experiments, reflections)
            refl_counts = flex.int()
            for expt_id in range(len(experiments)):
                refl_counts.append((refls_subset["id"] == expt_id).count(True))
            sort_order = flex.sort_permutation(refl_counts, reverse=True)
            if reflections.experiment_identifiers().keys():
                for idx in sorted(sort_order[: params.output.n_subset]):
                    subset_exp.append(experiments[idx])
                subset_refls = reflections.select(subset_exp)
                subset_refls.reset_ids()
            else:
                for expt_id, idx in enumerate(sort_order[: params.output.n_subset]):
                    subset_exp.append(experiments[idx])
                    refls = reflections.select(reflections["id"] == idx)
                    refls["id"] = flex.int(len(refls), expt_id)
                    subset_refls.extend(refls)

        experiments = subset_exp
        reflections = subset_refls

    # cluster the resulting experiments if requested
    if params.clustering.use and len(experiments) > 1:
        if params.clustering.max_clusters == 0:
            sys.exit("Error: max_clusters must be None or >0")
        clustered = do_unit_cell_clustering(
            experiments,
            reflections,
            dendrogram=params.clustering.dendrogram,
            threshold=params.clustering.threshold,
        )
        n_clusters = len(clustered)
        clusters = sorted(clustered.clusters, key=len, reverse=True)
        if params.clustering.max_clusters is not None:
            clusters = clusters[: params.clustering.max_clusters]
        if params.clustering.exclude_single_crystal_clusters:
            clusters = [c for c in clusters if len(c) > 1]
        clustered_experiments: List[ExperimentList] = [
            ExperimentList(
                [expt for expt in experiments if expt.identifier in cluster.lattice_ids]
            )
            for cluster in clusters
        ]
        clustered_reflections: List[flex.reflection_table] = [
            reflections.select_on_experiment_identifiers(cluster.lattice_ids)
            for cluster in clusters
        ]
        for i_cluster, (expts, refl) in enumerate(
            zip(clustered_experiments, clustered_reflections)
        ):
            refl.reset_ids()
            exp_name = os.path.splitext(params.output.experiments_filename)[0] + (
                "_cluster%d.expt" % (n_clusters - i_cluster)
            )
            refl_name = os.path.splitext(params.output.reflections_filename)[0] + (
                "_cluster%d.refl" % (n_clusters - i_cluster)
            )
            if params.output.max_batch_size is None:
                _save_output(
                    expts,
                    refl,
                    exp_name,
                    refl_name,
                )
            else:
                save_in_batches(
                    expts,
                    refl,
                    exp_name,
                    refl_name,
                    batch_size=params.output.max_batch_size,
                )

    else:
        if params.output.max_batch_size is None:
            _save_output(
                experiments,
                reflections,
                params.output.experiments_filename,
                params.output.reflections_filename,
            )
        else:
            save_in_batches(
                experiments,
                reflections,
                params.output.experiments_filename,
                params.output.reflections_filename,
                batch_size=params.output.max_batch_size,
            )
    return


def _save_output(experiments, reflections, exp_name, refl_name):
    # save output

    print(f"Saving combined experiments to {exp_name}")
    experiments.as_file(exp_name)
    print(f"Saving combined reflections to {refl_name}")
    reflections.as_file(refl_name)


def save_in_batches(experiments, reflections, exp_name, refl_name, batch_size=1000):
    for i, indices in enumerate(
        _split_equal_parts_of_length(
            list(range(len(experiments))), (len(experiments) // batch_size) + 1
        )
    ):
        batch_expts = ExperimentList()
        batch_refls = flex.reflection_table()
        if reflections.experiment_identifiers().keys():
            for sub_idx in indices:
                batch_expts.append(experiments[sub_idx])
            batch_refls = reflections.select(batch_expts)
            batch_refls.reset_ids()
        else:
            for sub_id, sub_idx in enumerate(indices):
                batch_expts.append(experiments[sub_idx])
                sub_refls = reflections.select(reflections["id"] == sub_idx)
                sub_refls["id"] = flex.int(len(sub_refls), sub_id)
                batch_refls.extend(sub_refls)
        exp_filename = os.path.splitext(exp_name)[0] + "_%03d.expt" % i
        ref_filename = os.path.splitext(refl_name)[0] + "_%03d.refl" % i
        _save_output(batch_expts, batch_refls, exp_filename, ref_filename)


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
        print("No Experiments found in the input")
        parser.print_help()
        return
    if not params.input.reflections:
        print("No reflection data found in the input")
        parser.print_help()
        return
    if len(params.input.reflections) != len(params.input.experiments):
        sys.exit(
            "The number of input reflections files does not match the "
            "number of input experiments"
        )

    try:
        flat_exps = flatten_experiments(params.input.experiments)
    except RuntimeError:
        sys.exit(
            "Unable to combine experiments. Are experiment IDs unique? "
            "You may need to run dials.assign_experiment_identifiers first to "
            "reset IDs."
        )

    run_with_preparsed(params, flat_exps)


if __name__ == "__main__":
    run()

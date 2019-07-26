#!/usr/bin/env dials.python
from __future__ import absolute_import, division, print_function

import os
import random

import dials.util
from dials.util import Sorry
from libtbx.phil import parse

help_message = """

Utility script to combine multiple reflections and experiments files into
one multi-experiment reflections and one experiments file. Experiments are
matched to reflections in the order they are provided as input.

Reference models can be chosen from any of the input experiments files. These
will replace all other models of that type in the output experiments file.
This is useful, for example, for combining mutiple experiments that should
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

    max_crystals = None
      .type = int
      .help = "Maximum number of crystals to cluster."

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
              "into seperate files, each at most max_batch_size number of"
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

    include scope dials.algorithms.integration.stills_significance_filter.phil_scope
  }
""",
    process_includes=True,
)


def find_experiment_in(experiment, all_experiments):
    """Search the phil experiment list and find where an experiment came from.

    :param Experiment experiment: The experiment to search for
    :param all_experiments:       The list of all experiments from phil
    :type  all_experiments:       list[dials.util.phil.FilenameDataWrapper[ExperimentList]]
    :returns:                     The filename and experiment ID
    :rtype:                       (str, int)
    """
    for source in all_experiments:
        try:
            experiment_list = list(source.data)
            index = experiment_list.index(experiment)
            return (source.filename, index)
        except ValueError:
            pass
    raise ValueError("Experiment not found")


class ComparisonError(Exception):
    """Exception to indicate problem with tolerance comparisons"""

    pass


class CombineWithReference(object):
    def __init__(
        self,
        beam=None,
        goniometer=None,
        scan=None,
        crystal=None,
        detector=None,
        params=None,
    ):

        self.ref_beam = beam
        self.ref_goniometer = goniometer
        self.ref_scan = scan
        self.ref_crystal = crystal
        self.ref_detector = detector
        self.tolerance = None
        self._last_imageset = None
        if params:
            if params.reference_from_experiment.compare_models:
                self.tolerance = params.reference_from_experiment.tolerance
            self.average_detector = params.reference_from_experiment.average_detector
        else:
            self.average_detector = False

    def __call__(self, experiment):
        from dxtbx.model.experiment_list import BeamComparison
        from dxtbx.model.experiment_list import DetectorComparison
        from dxtbx.model.experiment_list import GoniometerComparison
        from dxtbx.datablock import BeamDiff
        from dxtbx.datablock import DetectorDiff
        from dxtbx.datablock import GoniometerDiff

        if self.tolerance:
            compare_beam = BeamComparison(
                wavelength_tolerance=self.tolerance.beam.wavelength,
                direction_tolerance=self.tolerance.beam.direction,
                polarization_normal_tolerance=self.tolerance.beam.polarization_normal,
                polarization_fraction_tolerance=self.tolerance.beam.polarization_fraction,
            )
            compare_detector = DetectorComparison(
                fast_axis_tolerance=self.tolerance.detector.fast_axis,
                slow_axis_tolerance=self.tolerance.detector.slow_axis,
                origin_tolerance=self.tolerance.detector.origin,
            )
            compare_goniometer = GoniometerComparison(
                rotation_axis_tolerance=self.tolerance.goniometer.rotation_axis,
                fixed_rotation_tolerance=self.tolerance.goniometer.fixed_rotation,
                setting_rotation_tolerance=self.tolerance.goniometer.setting_rotation,
            )

        else:
            compare_beam = None
            compare_detector = None
            compare_goniometer = None

        if self.ref_beam:
            if compare_beam:
                if not compare_beam(self.ref_beam, experiment.beam):
                    diff = BeamDiff(
                        wavelength_tolerance=self.tolerance.beam.wavelength,
                        direction_tolerance=self.tolerance.beam.direction,
                        polarization_normal_tolerance=self.tolerance.beam.polarization_normal,
                        polarization_fraction_tolerance=self.tolerance.beam.polarization_fraction,
                    )
                    raise ComparisonError(
                        "\n".join(diff(self.ref_beam, experiment.beam))
                    )
            beam = self.ref_beam
        else:
            beam = experiment.beam

        if self.ref_detector and self.average_detector:
            detector = self.ref_detector
        elif self.ref_detector and not self.average_detector:
            if compare_detector:
                if not compare_detector(self.ref_detector, experiment.detector):
                    diff = DetectorDiff(
                        fast_axis_tolerance=self.tolerance.detector.fast_axis,
                        slow_axis_tolerance=self.tolerance.detector.slow_axis,
                        origin_tolerance=self.tolerance.detector.origin,
                    )
                    raise ComparisonError(
                        "\n".join(diff(self.ref_detector, experiment.detector))
                    )
            detector = self.ref_detector
        else:
            detector = experiment.detector

        if self.ref_goniometer:
            if compare_goniometer:
                if not compare_goniometer(self.ref_goniometer, experiment.goniometer):
                    diff = GoniometerDiff(
                        rotation_axis_tolerance=self.tolerance.goniometer.rotation_axis,
                        fixed_rotation_tolerance=self.tolerance.goniometer.fixed_rotation,
                        setting_rotation_tolerance=self.tolerance.goniometer.setting_rotation,
                    )
                    raise ComparisonError(
                        "\n".join(diff(self.ref_goniometer, experiment.goniometer))
                    )
            goniometer = self.ref_goniometer
        else:
            goniometer = experiment.goniometer

        if self.ref_scan:
            scan = self.ref_scan
        else:
            scan = experiment.scan

        if self.ref_crystal:
            crystal = self.ref_crystal
        else:
            crystal = experiment.crystal

        if self._last_imageset == experiment.imageset:
            imageset = self._last_imageset
        else:
            imageset = experiment.imageset
            self._last_imageset = imageset

        from dxtbx.model.experiment_list import Experiment

        return Experiment(
            beam=beam,
            detector=detector,
            scan=scan,
            goniometer=goniometer,
            crystal=crystal,
            imageset=imageset,
        )


class Cluster(object):
    def __init__(
        self, experiments, reflections, dendrogram=False, threshold=1000, n_max=None
    ):
        try:
            from xfel.clustering.cluster import Cluster
            from xfel.clustering.cluster_groups import unit_cell_info
        except ImportError:
            raise Sorry("clustering is not configured")
        import matplotlib.pyplot as plt

        ucs = Cluster.from_expts(
            refl_table=reflections, expts_list=experiments, n_images=n_max
        )
        self.clusters, _ = ucs.ab_cluster(
            threshold=threshold,
            log=True,  # log scale
            ax=plt.gca() if dendrogram else None,
            write_file_lists=False,
            schnell=False,
            doplot=dendrogram,
        )
        print(unit_cell_info(self.clusters))
        self.clustered_frames = {
            int(c.cname.split("_")[1]): c.members for c in self.clusters
        }
        if dendrogram:
            plt.tight_layout()
            plt.show()


class Script(object):
    def __init__(self):
        """Initialise the script."""
        from dials.util.options import OptionParser
        import libtbx.load_env

        # The script usage
        usage = (
            "usage: %s [options] [param.phil] "
            "experiments1.expt experiments2.expt reflections1.refl "
            "reflections2.refl..." % libtbx.env.dispatcher_name
        )

        # Create the parser
        self.parser = OptionParser(
            usage=usage,
            phil=phil_scope,
            read_reflections=True,
            read_experiments=True,
            check_format=False,
            epilog=help_message,
        )

    def run(self):
        """Execute the script."""
        params, options = self.parser.parse_args(show_diff_phil=True)
        self.run_with_preparsed(params, options)

    def run_with_preparsed(self, params, options):
        """Run combine_experiments, but allow passing in of parameters"""
        from dials.util.options import flatten_experiments

        # Try to load the models and data
        if len(params.input.experiments) == 0:
            print("No Experiments found in the input")
            self.parser.print_help()
            return
        if len(params.input.reflections) == 0:
            print("No reflection data found in the input")
            self.parser.print_help()
            return
        try:
            assert len(params.input.reflections) == len(params.input.experiments)
        except AssertionError:
            raise Sorry(
                "The number of input reflections files does not match the "
                "number of input experiments"
            )

        flat_exps = flatten_experiments(params.input.experiments)

        ref_beam = params.reference_from_experiment.beam
        ref_goniometer = params.reference_from_experiment.goniometer
        ref_scan = params.reference_from_experiment.scan
        ref_crystal = params.reference_from_experiment.crystal
        ref_detector = params.reference_from_experiment.detector

        if ref_beam is not None:
            try:
                ref_beam = flat_exps[ref_beam].beam
            except IndexError:
                raise Sorry("{} is not a valid experiment ID".format(ref_beam))

        if ref_goniometer is not None:
            try:
                ref_goniometer = flat_exps[ref_goniometer].goniometer
            except IndexError:
                raise Sorry("{} is not a valid experiment ID".format(ref_goniometer))

        if ref_scan is not None:
            try:
                ref_scan = flat_exps[ref_scan].scan
            except IndexError:
                raise Sorry("{} is not a valid experiment ID".format(ref_scan))

        if ref_crystal is not None:
            try:
                ref_crystal = flat_exps[ref_crystal].crystal
            except IndexError:
                raise Sorry("{} is not a valid experiment ID".format(ref_crystal))

        if ref_detector is not None:
            assert not params.reference_from_experiment.average_detector
            try:
                ref_detector = flat_exps[ref_detector].detector
            except IndexError:
                raise Sorry("{} is not a valid experiment ID".format(ref_detector))
        elif params.reference_from_experiment.average_detector:
            # Average all of the detectors together
            from scitbx.matrix import col

            def average_detectors(target, panelgroups, depth):
                # Recursive function to do the averaging

                if (
                    params.reference_from_experiment.average_hierarchy_level is None
                    or depth == params.reference_from_experiment.average_hierarchy_level
                ):
                    n = len(panelgroups)
                    sum_fast = col((0.0, 0.0, 0.0))
                    sum_slow = col((0.0, 0.0, 0.0))
                    sum_ori = col((0.0, 0.0, 0.0))

                    # Average the d matrix vectors
                    for pg in panelgroups:
                        sum_fast += col(pg.get_local_fast_axis())
                        sum_slow += col(pg.get_local_slow_axis())
                        sum_ori += col(pg.get_local_origin())
                    sum_fast /= n
                    sum_slow /= n
                    sum_ori /= n

                    # Re-orthagonalize the slow and the fast vectors by rotating around the cross product
                    c = sum_fast.cross(sum_slow)
                    a = sum_fast.angle(sum_slow, deg=True) / 2
                    sum_fast = sum_fast.rotate(c, a - 45, deg=True)
                    sum_slow = sum_slow.rotate(c, -(a - 45), deg=True)

                    target.set_local_frame(sum_fast, sum_slow, sum_ori)

                if target.is_group():
                    # Recurse
                    for i, target_pg in enumerate(target):
                        average_detectors(
                            target_pg, [pg[i] for pg in panelgroups], depth + 1
                        )

            ref_detector = flat_exps[0].detector
            average_detectors(
                ref_detector.hierarchy(), [e.detector.hierarchy() for e in flat_exps], 0
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
        from dials.array_family import flex

        reflections = flex.reflection_table()
        global_id = 0
        skipped_expts = 0
        from dxtbx.model.experiment_list import ExperimentList

        experiments = ExperimentList()

        # loop through the input, building up the global lists
        nrefs_per_exp = []
        for ref_wrapper, exp_wrapper in zip(
            params.input.reflections, params.input.experiments
        ):
            refs = ref_wrapper.data
            exps = exp_wrapper.data
            for i, exp in enumerate(exps):
                sel = refs["id"] == i
                sub_ref = refs.select(sel)
                n_sub_ref = len(sub_ref)
                if (
                    params.output.min_reflections_per_experiment is not None
                    and n_sub_ref < params.output.min_reflections_per_experiment
                ):
                    skipped_expts += 1
                    continue

                nrefs_per_exp.append(n_sub_ref)
                sub_ref["id"] = flex.int(len(sub_ref), global_id)
                if params.output.delete_shoeboxes and "shoebox" in sub_ref:
                    del sub_ref["shoebox"]
                reflections.extend(sub_ref)
                try:
                    experiments.append(combine(exp))
                except ComparisonError as e:
                    # When we failed tolerance checks, give a useful error message
                    (path, index) = find_experiment_in(exp, params.input.experiments)
                    raise Sorry(
                        "Model didn't match reference within required tolerance for experiment {} in {}:"
                        "\n{}\nAdjust tolerances or set compare_models=False to ignore differences.".format(
                            index, path, str(e)
                        )
                    )

                global_id += 1

        if (
            params.output.min_reflections_per_experiment is not None
            and skipped_expts > 0
        ):
            print(
                "Removed {0} experiments with fewer than {1} reflections".format(
                    skipped_expts, params.output.min_reflections_per_experiment
                )
            )

        # print number of reflections per experiment
        from libtbx.table_utils import simple_table

        header = ["Experiment", "Number of reflections"]
        rows = [(str(i), str(n)) for (i, n) in enumerate(nrefs_per_exp)]
        st = simple_table(rows, header)
        print(st.format())

        # save a random subset if requested
        if (
            params.output.n_subset is not None
            and len(experiments) > params.output.n_subset
        ):
            subset_exp = ExperimentList()
            subset_refls = flex.reflection_table()
            if params.output.n_subset_method == "random":
                n_picked = 0
                indices = list(range(len(experiments)))
                while n_picked < params.output.n_subset:
                    idx = indices.pop(random.randint(0, len(indices) - 1))
                    subset_exp.append(experiments[idx])
                    refls = reflections.select(reflections["id"] == idx)
                    refls["id"] = flex.int(len(refls), n_picked)
                    subset_refls.extend(refls)
                    n_picked += 1
                print(
                    "Selecting a random subset of {0} experiments out of {1} total.".format(
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
                    refl_counts.append(
                        len(refls_subset.select(refls_subset["id"] == expt_id))
                    )
                sort_order = flex.sort_permutation(refl_counts, reverse=True)
                for expt_id, idx in enumerate(sort_order[: params.output.n_subset]):
                    subset_exp.append(experiments[idx])
                    refls = reflections.select(reflections["id"] == idx)
                    refls["id"] = flex.int(len(refls), expt_id)
                    subset_refls.extend(refls)
                print(
                    "Selecting a subset of {0} experiments with highest number of reflections out of {1} total.".format(
                        params.output.n_subset, len(experiments)
                    )
                )

            elif params.output.n_subset_method == "significance_filter":
                from dials.algorithms.integration.stills_significance_filter import (
                    SignificanceFilter,
                )

                params.output.significance_filter.enable = True
                sig_filter = SignificanceFilter(params.output)
                refls_subset = sig_filter(experiments, reflections)
                refl_counts = flex.int()
                for expt_id in range(len(experiments)):
                    refl_counts.append(
                        len(refls_subset.select(refls_subset["id"] == expt_id))
                    )
                sort_order = flex.sort_permutation(refl_counts, reverse=True)
                for expt_id, idx in enumerate(sort_order[: params.output.n_subset]):
                    subset_exp.append(experiments[idx])
                    refls = reflections.select(reflections["id"] == idx)
                    refls["id"] = flex.int(len(refls), expt_id)
                    subset_refls.extend(refls)

            experiments = subset_exp
            reflections = subset_refls

        def save_in_batches(
            experiments, reflections, exp_name, refl_name, batch_size=1000
        ):
            from dxtbx.command_line.image_average import splitit

            for i, indices in enumerate(
                splitit(
                    list(range(len(experiments))), (len(experiments) // batch_size) + 1
                )
            ):
                batch_expts = ExperimentList()
                batch_refls = flex.reflection_table()
                for sub_id, sub_idx in enumerate(indices):
                    batch_expts.append(experiments[sub_idx])
                    sub_refls = reflections.select(reflections["id"] == sub_idx)
                    sub_refls["id"] = flex.int(len(sub_refls), sub_id)
                    batch_refls.extend(sub_refls)
                exp_filename = os.path.splitext(exp_name)[0] + "_%03d.expt" % i
                ref_filename = os.path.splitext(refl_name)[0] + "_%03d.refl" % i
                self._save_output(batch_expts, batch_refls, exp_filename, ref_filename)

        def combine_in_clusters(
            experiments_l, reflections_l, exp_name, refl_name, end_count
        ):
            result = []
            for cluster, experiment in enumerate(experiments_l):
                cluster_expts = ExperimentList()
                cluster_refls = flex.reflection_table()
                for i, expts in enumerate(experiment):
                    refls = reflections_l[cluster][i]
                    refls["id"] = flex.int(len(refls), i)
                    cluster_expts.append(expts)
                    cluster_refls.extend(refls)
                exp_filename = os.path.splitext(exp_name)[0] + (
                    "_cluster%d.expt" % (end_count - cluster)
                )
                ref_filename = os.path.splitext(refl_name)[0] + (
                    "_cluster%d.refl" % (end_count - cluster)
                )
                result.append(
                    (cluster_expts, cluster_refls, exp_filename, ref_filename)
                )
            return result

        # cluster the resulting experiments if requested
        if params.clustering.use:
            clustered = Cluster(
                experiments,
                reflections,
                dendrogram=params.clustering.dendrogram,
                threshold=params.clustering.threshold,
                n_max=params.clustering.max_crystals,
            )
            n_clusters = len(clustered.clustered_frames)

            def not_too_many(keeps):
                if params.clustering.max_clusters is not None:
                    return len(keeps) < params.clustering.max_clusters
                return True

            keep_frames = []
            sorted_keys = sorted(clustered.clustered_frames.keys())
            while len(clustered.clustered_frames) > 0 and not_too_many(keep_frames):
                keep_frames.append(clustered.clustered_frames.pop(sorted_keys.pop(-1)))
            if params.clustering.exclude_single_crystal_clusters:
                keep_frames = [k for k in keep_frames if len(k) > 1]
            clustered_experiments = [
                [f.experiment for f in frame_cluster] for frame_cluster in keep_frames
            ]
            clustered_reflections = [
                [f.reflections for f in frame_cluster] for frame_cluster in keep_frames
            ]
            list_of_combined = combine_in_clusters(
                clustered_experiments,
                clustered_reflections,
                params.output.experiments_filename,
                params.output.reflections_filename,
                n_clusters,
            )
            for saveable_tuple in list_of_combined:
                if params.output.max_batch_size is None:
                    self._save_output(*saveable_tuple)
                else:
                    save_in_batches(
                        *saveable_tuple, batch_size=params.output.max_batch_size
                    )
        else:
            if params.output.max_batch_size is None:
                self._save_output(
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

    def _save_output(self, experiments, reflections, exp_name, refl_name):
        # save output
        from dxtbx.model.experiment_list import ExperimentListDumper

        print("Saving combined experiments to {}".format(exp_name))
        dump = ExperimentListDumper(experiments)
        dump.as_json(exp_name)
        print("Saving combined reflections to {}".format(refl_name))
        reflections.as_pickle(refl_name)


if __name__ == "__main__":
    with dials.util.show_mail_on_error():
        script = Script()
        script.run()

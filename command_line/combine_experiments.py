from __future__ import absolute_import, division, print_function

import os
import random
import sys

import xfel.clustering.cluster
from dxtbx.command_line.image_average import splitit
from dxtbx.datablock import BeamDiff, DetectorDiff, GoniometerDiff
from dxtbx.model.experiment_list import (
    BeamComparison,
    DetectorComparison,
    Experiment,
    ExperimentList,
    GoniometerComparison,
)
from libtbx.phil import parse
from scitbx import matrix
from xfel.clustering.cluster_groups import unit_cell_info
from xfel.merging.application.input import file_lister

import dials.util
from dials.algorithms.integration.stills_significance_filter import SignificanceFilter
from dials.array_family import flex
from dials.util import tabulate
from dials.util.options import OptionParser, flatten_experiments
from dials.util.phil import ExperimentListConverters, ReflectionTableConverters

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

An MPI mode is available for multiprocessing of long file lists.

Examples::

  dials.combine_experiments experiments_0.expt experiments_1.expt \\
    reflections_0.refl reflections_1.refl \\
    reference_from_experiment.beam=0 \\
    reference_from_experiment.detector=0

  mpirun -n 32 dials.combine_experiments input.path=$PWD/out \\
    input.reflections_suffix=_strong.refl \\
    input.experiments_suffix=_imported.expt
"""

# The phil scope
phil_scope = parse(
    """
  input {
    path = None
      .type = str
      .multiple = True
      .help = A glob, directory, or file specifying experiments and reflections
      .help = to combine. Files name.expt and name.refl must be present as
      .help = matching pairs. This mode is only supported for MPI runs;
      .help = otherwise, supply expts and refls as command-line arguments.
    reflections_suffix = .refl
      .type = str
      .help = When input.path is given, find file names with this suffix for
      .help = reflections.
    experiments_suffix = .expt
      .type = str
      .help = When input.path is given, find file names with this suffix for
      .help = experiments..
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

        return Experiment(
            identifier=experiment.identifier,
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
        if dendrogram:
            import matplotlib.pyplot as plt

            axes = plt.gca()
        else:
            axes = None

        ucs = xfel.clustering.cluster.Cluster.from_expts(
            refl_table=reflections, expts_list=experiments, n_images=n_max
        )
        self.clusters, _ = ucs.ab_cluster(
            threshold=threshold,
            log=True,  # log scale
            ax=axes,
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
        # The script usage
        usage = (
            "usage: dials.combine_experiments [options] [param.phil] "
            "experiments1.expt experiments2.expt reflections1.refl "
            "reflections2.refl..."
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

    def run(self, args=None):
        """Execute the script."""

        try:
            from mpi4py import MPI

            comm = MPI.COMM_WORLD
            rank = comm.Get_rank()
            size = comm.Get_size()
        except ImportError:
            rank = 0
            size = 1

        if size == 1:
            params, options = self.parser.parse_args(args, show_diff_phil=True)
            assert not params.input.path, "input.path requires MPI support."
            self.run_with_preparsed(params, options)
            return

        if rank == 0:  # An MPI file loader!

            # Parse the command line
            params, options, fnames = self.parser.parse_args(
                args, show_diff_phil=False, quick_parse=True, return_unhandled=True
            )

            # Some basic validation
            i_argv = 0
            for fname in fnames:
                # Make sure the order of the arg list is unchanged
                try:
                    i_argv = sys.argv.index(fname, i_argv)
                except ValueError:
                    raise Exception("The order of the aruments was changed")

                # Make sure the only unhandled args are refls and expts
                if not (
                    fname.endswith(params.input.reflections_suffix)
                    or fname.endswith(params.input.experiments_suffix)
                ):
                    sys.exit("Did not understand argument %s" % fname)

            # Make (expt, refl) filename tuples
            if fnames:
                assert (
                    not params.input.path
                ), "Please provide path or filenames, not both"
                expts = [
                    f for f in fnames if f.endswith(params.input.experiments_suffix)
                ]
                refls = [
                    f for f in fnames if f.endswith(params.input.reflections_suffix)
                ]
                filenames = list(zip(expts, refls))
            else:
                lister = file_lister.file_lister(params)
                filenames = lister.filename_lister()

        if rank == 0:
            transmitted_info = params, options, fnames
        else:
            transmitted_info = None
        params, options, fnames = comm.bcast(transmitted_info, root=0)

        if rank == 0:  # Server

            results = [None] * len(filenames)
            working_ranks = [False] * size

            for i_fname, fname in enumerate(filenames):
                if i_fname % 1000 == 0:
                    print(i_fname)
                i_result, result, rankreq = comm.recv(source=MPI.ANY_SOURCE)
                if result:
                    results[i_result] = result
                comm.send((i_fname, fname), dest=rankreq)
                working_ranks[rankreq] = True

            for i in range(1, size):
                if not working_ranks[i]:
                    comm.send("done", dest=i)

            while working_ranks.count(True) > 0:
                i_result, result, rankreq = comm.recv(source=MPI.ANY_SOURCE)
                results[i_result] = result
                comm.send("done", dest=rankreq)
                working_ranks[rankreq] = False

        else:  # Client
            from libtbx.phil import strings_as_words as _words

            expt_converter = ExperimentListConverters(check_format=False)
            refl_converter = ReflectionTableConverters()
            result, i_result = None, None
            while True:
                comm.send((i_result, result, rank), dest=0)
                msg = comm.recv(source=0)
                if msg == "done":
                    break
                else:
                    i_fname, fname = msg
                    expt = expt_converter.from_words(_words([fname[0]]), None)
                    refl = refl_converter.from_words(_words([fname[1]]), None)
                    i_result = i_fname
                    result = (expt, refl)

        if rank != 0:
            exit()

        params.input.experiments = [r[0] for r in results if r[0]]
        params.input.reflections = [r[1] for r in results if r[1]]

        self.run_with_preparsed(params, options)

    def run_with_preparsed(self, params, options):
        """Run combine_experiments, but allow passing in of parameters"""
        # Try to load the models and data
        if not params.input.experiments:
            print("No Experiments found in the input")
            self.parser.print_help()
            return
        if not params.input.reflections:
            print("No reflection data found in the input")
            self.parser.print_help()
            return
        if len(params.input.reflections) != len(params.input.experiments):
            sys.exit(
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
                sys.exit("{} is not a valid experiment ID".format(ref_beam))

        if ref_goniometer is not None:
            try:
                ref_goniometer = flat_exps[ref_goniometer].goniometer
            except IndexError:
                sys.exit("{} is not a valid experiment ID".format(ref_goniometer))

        if ref_scan is not None:
            try:
                ref_scan = flat_exps[ref_scan].scan
            except IndexError:
                sys.exit("{} is not a valid experiment ID".format(ref_scan))

        if ref_crystal is not None:
            try:
                ref_crystal = flat_exps[ref_crystal].crystal
            except IndexError:
                sys.exit("{} is not a valid experiment ID".format(ref_crystal))

        if ref_detector is not None:
            assert not params.reference_from_experiment.average_detector
            try:
                ref_detector = flat_exps[ref_detector].detector
            except IndexError:
                sys.exit("{} is not a valid experiment ID".format(ref_detector))
        elif params.reference_from_experiment.average_detector:
            # Average all of the detectors together

            def average_detectors(target, panelgroups, depth):
                # Recursive function to do the averaging

                if (
                    params.reference_from_experiment.average_hierarchy_level is None
                    or depth == params.reference_from_experiment.average_hierarchy_level
                ):
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
                reflections.extend(sub_ref)
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

                global_id += 1

        if (
            params.output.min_reflections_per_experiment is not None
            and skipped_expts_min_refl > 0
        ):
            print(
                "Removed {0} experiments with fewer than {1} reflections".format(
                    skipped_expts_min_refl, params.output.min_reflections_per_experiment
                )
            )
        if (
            params.output.max_reflections_per_experiment is not None
            and skipped_expts_max_refl > 0
        ):
            print(
                "Removed {0} experiments with more than {1} reflections".format(
                    skipped_expts_max_refl, params.output.max_reflections_per_experiment
                )
            )

        # print number of reflections per experiment

        header = ["Experiment", "Number of reflections"]
        rows = [(str(i), str(n)) for (i, n) in enumerate(nrefs_per_exp)]
        print(tabulate(rows, header))

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
                    "Selecting a subset of {0} experiments with highest number of reflections out of {1} total.".format(
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

        def save_in_batches(
            experiments, reflections, exp_name, refl_name, batch_size=1000
        ):
            for i, indices in enumerate(
                splitit(
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
                    if refls.experiment_identifiers().keys():
                        identifier = refls.experiment_identifiers().values()[0]
                        id_val = refls.experiment_identifiers().keys()[0]
                        del refls.experiment_identifiers()[id_val]
                        refls["id"] = flex.int(len(refls), i)
                        refls.experiment_identifiers()[i] = identifier
                        refls.assert_experiment_identifiers_are_consistent(
                            experiment[i : i + 1]
                        )
                    else:
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
                ExperimentList([f.experiment for f in frame_cluster])
                for frame_cluster in keep_frames
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

        print("Saving combined experiments to {}".format(exp_name))
        experiments.as_file(exp_name)
        print("Saving combined reflections to {}".format(refl_name))
        reflections.as_file(refl_name)


@dials.util.show_mail_handle_errors()
def run(args=None):
    script = Script()
    script.run(args)


if __name__ == "__main__":
    run()

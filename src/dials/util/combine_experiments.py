from __future__ import annotations

import logging
import random
import sys
from dataclasses import dataclass

import dxtbx.model
import dxtbx.model.compare as compare
from dxtbx.model.experiment_list import (
    BeamComparison,
    DetectorComparison,
    Experiment,
    ExperimentList,
    GoniometerComparison,
)
from libtbx import phil
from scitbx import matrix

from dials.algorithms.clustering.unit_cell import cluster_unit_cells
from dials.algorithms.integration.stills_significance_filter import SignificanceFilter
from dials.algorithms.integration.stills_significance_filter import (
    phil_scope as sig_filter_phil_scope,
)
from dials.array_family import flex
from dials.util import tabulate

logger = logging.getLogger(__name__)


def find_experiment_in(
    experiment: Experiment, list_of_experiments: list[ExperimentList]
) -> tuple[int, int]:
    """Search the phil experiment list and find where an experiment came from."""
    for i, elist in enumerate(list_of_experiments):
        index = elist.find(experiment.identifier)
        if index == -1:
            continue
        return (i, index)
    raise ValueError("Experiment not found")


class ComparisonError(Exception):
    """Exception to indicate problem with tolerance comparisons"""

    pass


class CombineWithReference:
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
                    raise ComparisonError(
                        compare.beam_diff(
                            self.ref_beam,
                            experiment.beam,
                            wavelength_tolerance=self.tolerance.beam.wavelength,
                            direction_tolerance=self.tolerance.beam.direction,
                            polarization_normal_tolerance=self.tolerance.beam.polarization_normal,
                            polarization_fraction_tolerance=self.tolerance.beam.polarization_fraction,
                        )
                    )
            beam = self.ref_beam
        else:
            beam = experiment.beam

        if self.ref_detector and self.average_detector:
            detector = self.ref_detector
        elif self.ref_detector and not self.average_detector:
            if compare_detector:
                if not compare_detector(self.ref_detector, experiment.detector):
                    raise ComparisonError(
                        compare.detector_diff(
                            self.ref_detector,
                            experiment.detector,
                            fast_axis_tolerance=self.tolerance.detector.fast_axis,
                            slow_axis_tolerance=self.tolerance.detector.slow_axis,
                            origin_tolerance=self.tolerance.detector.origin,
                        )
                    )
            detector = self.ref_detector
        else:
            detector = experiment.detector

        if self.ref_goniometer:
            if compare_goniometer:
                if not compare_goniometer(self.ref_goniometer, experiment.goniometer):
                    raise ComparisonError(
                        compare.goniometer_diff(
                            self.ref_goniometer,
                            experiment.goniometer,
                            rotation_axis_tolerance=self.tolerance.goniometer.rotation_axis,
                            fixed_rotation_tolerance=self.tolerance.goniometer.fixed_rotation,
                            setting_rotation_tolerance=self.tolerance.goniometer.setting_rotation,
                        )
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


def do_unit_cell_clustering(experiments, dendrogram=False, threshold=1000):
    if dendrogram:
        import matplotlib.pyplot as plt

        ax = plt.gca()
    else:
        ax = None

    crystal_symmetries = [expt.crystal.get_crystal_symmetry() for expt in experiments]
    clustering = cluster_unit_cells(
        crystal_symmetries,
        lattice_ids=list(experiments.identifiers()),
        threshold=threshold,
        ax=ax,
        no_plot=not dendrogram,
    )

    if dendrogram:
        ax.set_yscale("symlog", linthresh=1)
        plt.tight_layout()
        plt.show()

    return clustering


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
    flat_exps: ExperimentList, reference_from_experiment: phil.scope_extract
) -> tuple[
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


@dataclass
class significance_params:
    significance_filter: phil.scope_extract


def _select_random_subset(experiments, reflections=None, n_subset: int = 1):
    n_picked = 0
    indices = list(range(len(experiments)))
    subset_exp = ExperimentList()
    subset_refls = flex.reflection_table()
    if reflections:
        if reflections.experiment_identifiers().keys():
            indices_to_sel = []
            while n_picked < n_subset:
                idx = indices.pop(random.randint(0, len(indices) - 1))
                indices_to_sel.append(idx)
                n_picked += 1
            # make sure select in order.
            for idx in sorted(indices_to_sel):
                subset_exp.append(experiments[idx])
            subset_refls = reflections.select(subset_exp)
            subset_refls.reset_ids()
        else:
            while n_picked < n_subset:
                idx = indices.pop(random.randint(0, len(indices) - 1))
                subset_exp.append(experiments[idx])
                refls = reflections.select(reflections["id"] == idx)
                refls["id"] = flex.int(len(refls), n_picked)
                subset_refls.extend(refls)
                n_picked += 1
        logger.info(
            f"Selecting a random subset of {n_subset} experiments out of {len(experiments)} total."
        )
    else:
        while n_picked < n_subset:
            idx = indices.pop(random.randint(0, len(indices) - 1))
            subset_exp.append(experiments[idx])
            n_picked += 1
        logger.info(
            f"Selecting a random subset of {n_subset} experiments out of {len(experiments)} total."
        )
    return subset_exp, subset_refls


def select_subset(
    experiments: ExperimentList,
    reflections: flex.reflection_table,
    n_subset: int = 1,
    n_subset_method: str = "random",
    n_refl_panel_list=None,
    significance_filter: phil.scope_extract | None = None,
) -> tuple[ExperimentList, flex.reflection_table]:
    subset_exp = ExperimentList()
    subset_refls = flex.reflection_table()
    if n_subset_method == "random":  # Doesn't require reflections
        subset_exp, subset_refls = _select_random_subset(
            experiments, reflections, n_subset
        )
    elif n_subset_method == "n_refl":
        if n_refl_panel_list is None:
            refls_subset = reflections
        else:
            sel = flex.bool(len(reflections), False)
            for p in n_refl_panel_list:
                sel |= reflections["panel"] == p
            refls_subset = reflections.select(sel)
        refl_counts = flex.int()
        for expt_id in range(len(experiments)):
            refl_counts.append((refls_subset["id"] == expt_id).count(True))
        sort_order = flex.sort_permutation(refl_counts, reverse=True)
        if reflections.experiment_identifiers().keys():
            for idx in sorted(sort_order[:n_subset]):
                subset_exp.append(experiments[idx])
            subset_refls = reflections.select(subset_exp)
            subset_refls.reset_ids()
        else:
            for expt_id, idx in enumerate(sort_order[:n_subset]):
                subset_exp.append(experiments[idx])
                refls = reflections.select(reflections["id"] == idx)
                refls["id"] = flex.int(len(refls), expt_id)
                subset_refls.extend(refls)
        logger.info(
            f"Selecting a subset of {n_subset} experiments with highest number of reflections out of {len(experiments)} total."
        )

    elif n_subset_method == "significance_filter":
        if significance_filter is None:
            significance_filter = sig_filter_phil_scope.extract()
        significance_filter.enable = True
        sig_filter = SignificanceFilter(significance_params(significance_filter))
        refls_subset = sig_filter(experiments, reflections)
        refl_counts = flex.int()
        for expt_id in range(len(experiments)):
            refl_counts.append((refls_subset["id"] == expt_id).count(True))
        sort_order = flex.sort_permutation(refl_counts, reverse=True)
        if reflections.experiment_identifiers().keys():
            for idx in sorted(sort_order[:n_subset]):
                subset_exp.append(experiments[idx])
            subset_refls = reflections.select(subset_exp)
            subset_refls.reset_ids()
        else:
            for expt_id, idx in enumerate(sort_order[:n_subset]):
                subset_exp.append(experiments[idx])
                refls = reflections.select(reflections["id"] == idx)
                refls["id"] = flex.int(len(refls), expt_id)
                subset_refls.extend(refls)
    else:
        raise ValueError(f"Invalid option for n_subset_method: {n_subset_method}")

    return subset_exp, subset_refls


def combine_experiments_no_reflections(params, experiment_lists):
    """Run combine_experiments, without corresponding reflection tables"""
    flat_exps = ExperimentList()
    try:
        for elist in experiment_lists:
            flat_exps.extend(elist)
    except RuntimeError:
        sys.exit(  # FIXME
            "Unable to combine experiments. Are experiment IDs unique? "
            "You may need to run dials.assign_experiment_identifiers first to "
            "reset IDs."
        )

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
    experiments = ExperimentList()
    for elist in experiment_lists:
        for expt in elist:
            try:
                experiments.append(combine(expt))
            except ComparisonError as e:
                # When we failed tolerance checks, give a useful error message
                (i, index) = find_experiment_in(expt, experiment_lists)
                sys.exit(  # FIXME - raise RuntimeError?
                    f"Model didn't match reference within required tolerance for experiment {index} in input file {i}:"
                    f"\n{str(e)}\nAdjust tolerances or set compare_models=False to ignore differences."
                )
    # select a subset if requested
    if params.output.n_subset is not None and len(experiments) > params.output.n_subset:
        assert params.output.n_subset_method == "random", (
            "Combining only experiments and not reflections is only possible with n_subset_method=random if n_subset is not None"
        )
        experiments, _ = _select_random_subset(
            experiments,
            reflections=None,
            n_subset=params.output.n_subset,
        )
    return experiments


def combine_experiments(params, experiment_lists, reflection_tables):
    """Run combine_experiments"""

    flat_exps = ExperimentList()
    try:
        for elist in experiment_lists:
            flat_exps.extend(elist)
    except RuntimeError:
        sys.exit(  # FIXME
            "Unable to combine experiments. Are experiment IDs unique? "
            "You may need to run dials.assign_experiment_identifiers first to "
            "reset IDs."
        )

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
    experiments = ExperimentList()
    global_id = 0
    skipped_expts_min_refl = 0
    skipped_expts_max_refl = 0

    # loop through the input, building up the global lists
    nrefs_per_exp = []
    for refs, exps in zip(reflection_tables, experiment_lists):
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
                (i, index) = find_experiment_in(exp, experiment_lists)
                sys.exit(  # FIXME - raise RuntimeError?
                    f"Model didn't match reference within required tolerance for experiment {index} in input file {i}:"
                    f"\n{str(e)}\nAdjust tolerances or set compare_models=False to ignore differences."
                )

            # Rewrite imageset_id, if the experiment has an imageset
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

    # Finished building global lists

    if (
        params.output.min_reflections_per_experiment is not None
        and skipped_expts_min_refl > 0
    ):
        logger.info(
            f"Removed {skipped_expts_min_refl} experiments with fewer than {params.output.min_reflections_per_experiment} reflections"
        )
    if (
        params.output.max_reflections_per_experiment is not None
        and skipped_expts_max_refl > 0
    ):
        logger.info(
            f"Removed {skipped_expts_max_refl} experiments with more than {params.output.max_reflections_per_experiment} reflections"
        )

    # print number of reflections per experiment

    header = ["Experiment", "Number of reflections"]
    rows = [(str(i), str(n)) for (i, n) in enumerate(nrefs_per_exp)]
    logger.info(tabulate(rows, header))

    # select a subset if requested
    if params.output.n_subset is not None and len(experiments) > params.output.n_subset:
        experiments, reflections = select_subset(
            experiments,
            reflections,
            n_subset=params.output.n_subset,
            n_subset_method=params.output.n_subset_method,
            n_refl_panel_list=params.output.n_refl_panel_list,
            significance_filter=params.output.significance_filter,
        )

    return experiments, reflections

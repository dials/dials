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
            # For XFEL stills, per-frame wavelengths live in the scan
            # "wavelength" property — they are not a beam attribute. After
            # decode()'s XFELBeam → monochromatic Beam auto-resolve, each
            # experiment's beam carries its own frame wavelength, so a
            # strict wavelength_tolerance would always fail. Relax it to
            # infinity when both reference and current experiment are
            # stills; geometry (direction, polarization) is still checked.
            wavelength_tolerance = self.tolerance.beam.wavelength
            ref_is_still = self.ref_scan is not None and self.ref_scan.is_still()
            cur_is_still = experiment.scan is not None and experiment.scan.is_still()
            if ref_is_still and cur_is_still:
                wavelength_tolerance = float("inf")
            compare_beam = BeamComparison(
                wavelength_tolerance=wavelength_tolerance,
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
            # For stills, each experiment carries its own wavelength in the scan
            # properties. Replacing with the reference scan would clobber those
            # per-frame wavelengths. Keep the experiment's own scan; serialisation
            # automatically consolidates stills scans on write.
            if (
                self.ref_scan.is_still()
                and experiment.scan is not None
                and experiment.scan.is_still()
            ):
                scan = experiment.scan
            else:
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
            profile=experiment.profile,
            imageset=imageset,
            scaling_model=experiment.scaling_model,
            history=experiment.history,
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


def consolidate_stills_imagesets(experiments, reflections=None):
    """Merge experiments from the same source file into one shared ImageSequence.

    When N per-worker .expt files are combined, each file contributes its own
    ImageSequence Python object even when all files reference the same source
    data file. This causes to_dict() to write N imageset entries. Group by
    source path and re-create one shared ImageSequence per unique path.

    Also prunes a shared ImageSequence back to its surviving experiments after a
    tool (combine filtering, split_still_data) drops some experiments, so the
    consolidated scan's single_file_indices does not over-count the scan_points
    on read.

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
    # Skip rotation data. Experiment.is_still() is the full stills contract
    # (scan-None / scan.is_still() / goniometer-None, and correctly excludes
    # Laue and ToF), not just the imageset scan's oscillation width.
    if not experiments[0].is_still():
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

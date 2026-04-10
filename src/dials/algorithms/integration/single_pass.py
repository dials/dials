"""Phase 3 MVP: single-pass chunked modeller-integrator skeleton."""

from __future__ import annotations

import logging

from dials.algorithms.shoebox import MaskCode
from dials_algorithms_integration_integrator_ext import Executor

logger = logging.getLogger(__name__)


class SliceReadinessTracker:
    """Tracks finalized profile-model z-slices.

    A cell is ready for fitting once its slice has been finalized via
    mark_slice_finalized().  Zero cross-cell coupling: finalizing slice k
    leaves other slices in their pre-finalized (additive) state.
    """

    def __init__(self, num_slices: int, cells_per_slice: int) -> None:
        self.num_slices = num_slices
        self.cells_per_slice = cells_per_slice
        self.finalized_slices: set[int] = set()

    def is_cell_ready(self, cell_idx: int) -> bool:
        return (cell_idx // self.cells_per_slice) in self.finalized_slices

    def mark_slice_finalized(self, slice_idx: int) -> None:
        self.finalized_slices.add(slice_idx)


class PendingFitQueue:
    """Holds reflections whose profile-model cells are not yet finalized.

    Keyed by nearest-cell slice index.  Each entry stores a deep copy of the
    sub-table row (with shoebox intact) and the master row index.  The deep
    copy is essential: ShoeboxProcessor deallocates shoeboxes after
    process() returns, so we must capture them inside _handle_fit() before
    returning.
    """

    def __init__(self, cells_per_slice: int = 9) -> None:
        self._by_slice: dict[int, list] = {}  # slice_idx -> [(master_idx, 1-row table)]
        self._cells_per_slice = cells_per_slice

    def enqueue(self, pending_subtable, nearest_cells, cells_per_slice: int) -> None:
        """Store pending reflections keyed by their nearest cell's slice index.

        pending_subtable is a deep-copied reflection_table with shoeboxes.
        nearest_cells is a list/array of cell indices (one per row).
        """
        master_indices = list(pending_subtable["_master_idx"])
        for i in range(len(pending_subtable)):
            cell = nearest_cells[i]
            slice_idx = cell // cells_per_slice
            if slice_idx not in self._by_slice:
                self._by_slice[slice_idx] = []
            # Store the master index and a 1-row sub-table
            row = pending_subtable[i : i + 1]
            self._by_slice[slice_idx].append((master_indices[i], row))

    def drain_by_slice(self, slice_idx: int, profile_fitter, master_table) -> None:
        """Fit all pending reflections for this slice and write results to master.

        profile_fitter must be the raw MultiExpProfileModeller (not the
        ValidatedMultiExpProfileModeller): compute_fitted_intensity calls
        .fit() which exists on MultiExpProfileModeller but NOT on
        ValidatedMultiExpProfileModeller.
        """
        from dials.array_family import flex

        entries = self._by_slice.pop(slice_idx, [])
        if not entries:
            return

        # Concatenate all pending rows for this slice into one table
        combined = entries[0][1]
        master_idxs = [entries[0][0]]
        for master_idx, row_table in entries[1:]:
            combined.extend(row_table)
            master_idxs.append(master_idx)

        # Fit the combined table (profile_fitter is the raw MultiExpProfileModeller)
        combined.compute_fitted_intensity(profile_fitter)

        # Write fit results back to master table
        idx_array = flex.size_t(master_idxs)
        for col in [
            "intensity.prf.value",
            "intensity.prf.variance",
            "profile.correlation",
        ]:
            if col in combined:
                master_table[col].set_selected(idx_array, combined[col])

        # Propagate IntegratedPrf flag
        prf_flag = combined.get_flags(combined.flags.integrated_prf)
        for i in range(len(idx_array)):
            if prf_flag[i]:
                master_table["flags"][idx_array[i]] |= combined.flags.integrated_prf

    def flush_all(self, profile_fitter, master_table) -> None:
        """Drain all remaining entries (for tail slices)."""
        for slice_idx in sorted(self._by_slice.keys()):
            self.drain_by_slice(slice_idx, profile_fitter, master_table)

    def is_empty(self) -> bool:
        return len(self._by_slice) == 0


class SinglePassExecutor(Executor):
    """Combines profile modelling + fitting in one process() call.

    Step order replicates IntegratorExecutor.process() lines 897-935; profile
    modelling is inserted after summed-intensity. Steps 1-8 run on the full batch;
    step 9 extracts the reference subset for modelling only (Q8/integrator.py:1013).
    """

    __getstate_manages_dict__ = 1

    def __init__(
        self,
        experiments,
        profile_fitter,
        sampler,
        readiness_tracker: SliceReadinessTracker,
        pending_queue: PendingFitQueue,
        valid_foreground_threshold: float = 0.75,
    ) -> None:
        self.experiments = experiments
        self.profile_fitter = profile_fitter
        self.sampler = sampler
        self.readiness_tracker = readiness_tracker
        self.pending_queue = pending_queue
        self.valid_foreground_threshold = valid_foreground_threshold
        super().__init__()

    def initialize(self, frame0: int, frame1: int, reflections) -> None:
        """Mirrors IntegratorExecutor.initialize(): find overlaps + log."""
        reflections.find_overlaps(self.experiments)
        logger.info(
            " Single-pass chunk: frames %d -> %d, %d reflections",
            frame0,
            frame1,
            len(reflections),
        )

    def process(self, frame: int, reflections) -> None:
        """Process one frame batch: modelling + fitting in one pass."""
        # Steps 1-3: overload / shoebox mask / invalid-pixel flags
        reflections.is_overloaded(self.experiments)
        reflections.compute_mask(self.experiments)
        reflections.contains_invalid_pixels()

        # Step 4: fraction_valid gate -> dont_integrate (lines 906-913)
        sbox = reflections["shoebox"]
        nvalfg = sbox.count_mask_values(MaskCode.Valid | MaskCode.Foreground)
        nforeg = sbox.count_mask_values(MaskCode.Foreground)
        reflections.set_flags(
            nvalfg.as_double() / nforeg.as_double() < self.valid_foreground_threshold,
            reflections.flags.dont_integrate,
        )

        # Steps 5-7: background, centroid, summed intensity (full batch)
        reflections.compute_background(self.experiments)
        reflections.compute_centroid(self.experiments)
        reflections.compute_summed_intensity()

        # Unset integrated_sum for dont_integrate reflections.  The two-pass
        # path's ReflectionLookup (integrator.h:400) excludes these from all
        # job assignments so they never reach compute_summed_intensity at all.
        di_mask = reflections.get_flags(reflections.flags.dont_integrate)
        reflections.unset_flags(di_mask, reflections.flags.integrated_sum)

        # Step 8: num_pixels.* accounting (lines 928-935)
        reflections["num_pixels.valid"] = sbox.count_mask_values(MaskCode.Valid)
        reflections["num_pixels.background"] = sbox.count_mask_values(
            MaskCode.Valid | MaskCode.Background
        )
        reflections["num_pixels.background_used"] = sbox.count_mask_values(
            MaskCode.Valid | MaskCode.Background | MaskCode.BackgroundUsed
        )
        reflections["num_pixels.foreground"] = nvalfg

        # Steps 9-10: reference subset (row-order-preserving) -> modelling only.
        # Guard: production MultiExpProfileModeller.model() asserts size > 0.
        ref_subset = reflections.select(
            reflections.get_flags(reflections.flags.reference_spot)
        )
        if len(ref_subset) > 0:
            self.profile_fitter.model(ref_subset)

        # Step 11: route reflections to inline-fit or pending queue
        self._handle_fit(frame, reflections)

        nsum = reflections.get_flags(reflections.flags.integrated_sum).count(True)
        nprf = reflections.get_flags(reflections.flags.integrated_prf).count(True)
        logger.debug(
            " Single-pass frame %d: %d sum + %d prf / %d",
            frame + 1,
            nsum,
            nprf,
            len(reflections),
        )

    def _handle_fit(self, frame: int, reflections) -> None:
        """Split reflections into inline-fit (cell ready) vs pending (cell not ready).

        For ready reflections: fit immediately via compute_fitted_intensity.
        For pending reflections: deep-copy with shoeboxes and enqueue.
        Must complete before returning — ShoeboxProcessor deallocates shoeboxes after.

        OPEN-F guard: reflections with z outside scan_range cause GridSampler.nearest()
        to assert-fail.  These 19 boundary reflections are skipped (they are also
        filtered by check1 in production and never reach fitting there).
        """
        from dials.array_family import flex

        scan_range = self.experiments[0].scan.get_array_range()
        xyz = reflections["xyzcal.px"]
        ready_mask = flex.bool(len(reflections), False)
        nearest_cells = flex.int(len(reflections), -1)  # -1 = boundary/error

        for i in range(len(reflections)):
            z = xyz[i][2]
            # Guard: GridSampler.nearest() asserts z in [scan_lo, scan_hi)
            if z < scan_range[0] or z >= scan_range[1]:
                continue  # boundary reflection — skip routing, let it fall unfitted
            cell = self.sampler.nearest(0, xyz[i])
            nearest_cells[i] = cell
            if self.readiness_tracker.is_cell_ready(cell):
                ready_mask[i] = True

        # Inline-fit the ready subset
        if ready_mask.count(True) > 0:
            ready_subset = reflections.select(ready_mask)
            ready_subset.compute_fitted_intensity(self.profile_fitter.modellers[0])
            # Write fit results back to the sub-table; ShoeboxProcessor propagates
            # these to master via set_selected_rows_index after process() returns.
            reflections.set_selected(ready_mask, ready_subset)

        # Enqueue the pending subset — deep copy preserves shoebox data before
        # ShoeboxProcessor deallocates them on return from process().
        pending_mask = ~ready_mask & (nearest_cells >= 0)
        if pending_mask.count(True) > 0:
            pending_subset = reflections.select(pending_mask).copy()  # deep copy!
            pending_cells = nearest_cells.select(pending_mask)
            self.pending_queue.enqueue(
                pending_subset,
                list(pending_cells),
                self.readiness_tracker.cells_per_slice,
            )

    def finalize(self) -> None:
        pass

    def data(self):
        return self.profile_fitter


class SinglePassBlockExecutor(Executor):
    """Executor for block-parallel single-pass integration (nproc>1).

    Unlike SinglePassExecutor, this class does NOT route reflections to
    inline-fit or pending queue during process().  It only performs steps 1-9
    (mask, bg, centroid, sumint, model) and defers all profile fitting to after
    the frame loop (handled by SinglePassBlockTask).

    No sampler, readiness_tracker, or pending_queue dependencies.
    """

    __getstate_manages_dict__ = 1

    def __init__(
        self,
        experiments,
        profile_fitter,
        valid_foreground_threshold: float = 0.75,
    ) -> None:
        self.experiments = experiments
        self.profile_fitter = profile_fitter
        self.valid_foreground_threshold = valid_foreground_threshold
        super().__init__()

    def initialize(self, frame0: int, frame1: int, reflections) -> None:
        """Mirrors SinglePassExecutor.initialize(): find overlaps + log."""
        reflections.find_overlaps(self.experiments)
        logger.info(
            " Block executor: frames %d -> %d, %d reflections",
            frame0,
            frame1,
            len(reflections),
        )

    def process(self, frame: int, reflections) -> None:
        """Process one frame batch: steps 1-9 only (no profile fitting).

        Identical to SinglePassExecutor.process() lines 152-192 minus the
        _handle_fit() call.  Profile fitting is deferred to after the frame
        loop in SinglePassBlockTask.__call__().
        """
        # Steps 1-3: overload / shoebox mask / invalid-pixel flags
        reflections.is_overloaded(self.experiments)
        reflections.compute_mask(self.experiments)
        reflections.contains_invalid_pixels()

        # Step 4: fraction_valid gate -> dont_integrate
        sbox = reflections["shoebox"]
        nvalfg = sbox.count_mask_values(MaskCode.Valid | MaskCode.Foreground)
        nforeg = sbox.count_mask_values(MaskCode.Foreground)
        reflections.set_flags(
            nvalfg.as_double() / nforeg.as_double() < self.valid_foreground_threshold,
            reflections.flags.dont_integrate,
        )

        # Steps 5-7: background, centroid, summed intensity (full batch)
        reflections.compute_background(self.experiments)
        reflections.compute_centroid(self.experiments)
        reflections.compute_summed_intensity()

        # Unset integrated_sum for dont_integrate reflections.
        di_mask = reflections.get_flags(reflections.flags.dont_integrate)
        reflections.unset_flags(di_mask, reflections.flags.integrated_sum)

        # Step 8: num_pixels.* accounting
        reflections["num_pixels.valid"] = sbox.count_mask_values(MaskCode.Valid)
        reflections["num_pixels.background"] = sbox.count_mask_values(
            MaskCode.Valid | MaskCode.Background
        )
        reflections["num_pixels.background_used"] = sbox.count_mask_values(
            MaskCode.Valid | MaskCode.Background | MaskCode.BackgroundUsed
        )
        reflections["num_pixels.foreground"] = nvalfg

        # Step 9: reference subset -> modelling only (no fitting here)
        ref_subset = reflections.select(
            reflections.get_flags(reflections.flags.reference_spot)
        )
        if len(ref_subset) > 0:
            self.profile_fitter.model(ref_subset)

    def finalize(self) -> None:
        pass

    def data(self):
        return self.profile_fitter

    def __getinitargs__(self):
        return (self.experiments, self.profile_fitter, self.valid_foreground_threshold)


class SinglePassBlockTask:
    """Per-worker task for nproc>1 single-pass integration.

    Builds a local profile model, runs the frame loop over the extended range,
    finalizes all model cells, fits all reflections, then filters to the core
    range.  Returns a dials.algorithms.integration.Result.
    """

    def __init__(
        self,
        index,
        experiments,
        reflections,
        extended_range,
        core_range,
        num_scan_points,
    ):
        self.index = index
        self.experiments = experiments
        self.reflections = reflections
        self.extended_range = extended_range  # (frame0, frame1)
        self.core_range = core_range  # (frame0, frame1)
        self.num_scan_points = num_scan_points

    def __call__(self):
        from time import time

        import dials.algorithms.integration
        from dials.algorithms.integration.processor import ShoeboxProcessor
        from dials.algorithms.integration.validation import (
            ValidatedMultiExpProfileModeller,
        )
        from dials.algorithms.profile_model.modeller import MultiExpProfileModeller
        from dials.array_family import flex
        from dials.model.data import make_image

        start_time = time()

        # Build local profile modeller (same pattern as ChunkDriver._build_validated_modeller)
        pf = MultiExpProfileModeller()
        for expr in self.experiments:
            pf.add(_make_modeller(expr, self.num_scan_points))
        profile_fitter = ValidatedMultiExpProfileModeller([pf])

        # Partiality is computed by run_single_pass_parallel() before dispatch.

        # Pre-create fitting output columns (needed for sentinel check below;
        # compute_fitted_intensity creates them but leaves unfitted reflections
        # with undefined values)
        n = len(self.reflections)
        for col in (
            "intensity.prf.value",
            "intensity.prf.variance",
            "profile.correlation",
        ):
            if col not in self.reflections:
                self.reflections[col] = flex.double(n, 0.0)

        # Create executor
        executor = SinglePassBlockExecutor(self.experiments, profile_fitter)

        # Slice imageset to extended range (same pattern as processor.Task.__call__)
        imageset = self.experiments[0].imageset
        allowed_range = imageset.get_array_range()
        index0 = self.extended_range[0] - allowed_range[0]
        index1 = self.extended_range[1] - allowed_range[0]
        imageset = imageset[index0:index1]
        frame0, frame1 = imageset.get_array_range()

        # Initialize executor
        executor.initialize(frame0, frame1, self.reflections)

        # Allocate shoeboxes (not allocated yet)
        self.reflections["shoebox"] = flex.shoebox(
            self.reflections["panel"],
            self.reflections["bbox"],
            allocate=False,
            flatten=False,
        )

        # Create ShoeboxProcessor with save=True to keep shoeboxes for fitting
        n_panels = len(imageset.get_detector())
        proc = ShoeboxProcessor(self.reflections, n_panels, frame0, frame1, True)

        # Frame loop
        import sys

        n_frames = len(imageset)
        log_interval = max(1, n_frames // 10)
        for i in range(n_frames):
            image = imageset.get_corrected_data(i)
            mask = imageset.get_mask(i)
            proc.next(make_image(image, mask), executor)
            del image, mask
            if (i + 1) % log_interval == 0 or i + 1 == n_frames:
                print(
                    f" Worker {self.index}: {i + 1}/{n_frames} frames",
                    file=sys.stderr,
                    flush=True,
                )
        assert proc.finished(), "ShoeboxProcessor did not finish"

        # Finalize populated model cells (empty cells have uninitialized data
        # arrays and will assert-fail in finalize_cell)
        modeller = profile_fitter.modellers[0][0]  # GaussianRSProfileModeller
        for cell_idx in range(modeller.size()):
            if modeller.n_reflections(cell_idx) > 0:
                modeller.finalize_cell(cell_idx)

        # Guard: reflections with xyzcal.px.z outside scan range crash
        # GridSampler.nearest() (hard assert). Mark them dont_integrate so
        # compute_fitted_intensity skips them; the sentinel fixup below gives
        # them variance=-1 and failed_during_profile_fitting, matching two-pass.
        scan_range = self.experiments[0].scan.get_array_range()
        z_vals = self.reflections["xyzcal.px"].parts()[2]
        oob = (z_vals < scan_range[0]) | (z_vals >= scan_range[1])
        self.reflections.set_flags(oob, self.reflections.flags.dont_integrate)

        # Fit ALL reflections using the raw MultiExpProfileModeller
        self.reflections.compute_fitted_intensity(profile_fitter.modellers[0])

        # Clear the temporary dont_integrate guard so it doesn't pollute flags
        self.reflections.unset_flags(oob, self.reflections.flags.dont_integrate)

        # Deallocate shoeboxes
        del self.reflections["shoebox"]

        # Boundary sentinel fixup (same as ChunkDriver.run() lines 423-429)
        prf_var = self.reflections["intensity.prf.variance"]
        prf_flag = self.reflections.get_flags(self.reflections.flags.integrated_prf)
        never_fitted = (prf_var == 0.0) & ~prf_flag
        prf_var.set_selected(never_fitted, -1.0)
        self.reflections.set_flags(
            never_fitted, self.reflections.flags.failed_during_profile_fitting
        )

        # Filter to core range using xyzcal.px z coordinate.
        # First worker also claims reflections with z < scan_lo (scan-start
        # partials whose predicted centroid is before frame 0).
        xyz = self.reflections["xyzcal.px"]
        z_vals = xyz.parts()[2]
        core_lo, core_hi = self.core_range
        orphan = (
            z_vals < scan_range[0]
            if core_lo == scan_range[0]
            else flex.bool(len(z_vals), False)
        )
        if core_hi == scan_range[1]:
            keep = (z_vals >= core_lo) | orphan
        else:
            keep = ((z_vals >= core_lo) & (z_vals < core_hi)) | orphan
        self.reflections = self.reflections.select(keep)

        # Finalize executor
        executor.finalize()

        return dials.algorithms.integration.Result(
            index=self.index,
            reflections=self.reflections,
            data=profile_fitter,
            read_time=0,
            extract_time=proc.extract_time(),
            process_time=proc.process_time(),
            total_time=time() - start_time,
        )


def _build_sampler(experiments, num_scan_points: int):
    """Construct GridSampler matching production defaults (scan_step=5)."""
    from dials.algorithms.profile_model.modeller import GridSampler

    return GridSampler(
        experiments[0].detector[0].get_image_size(),
        experiments[0].scan.get_array_range(),
        (3, 3, num_scan_points),
    )


def _make_modeller(expr, num_scan_points: int):
    """Construct GaussianRSProfileModeller with production defaults.
    expr.profile.params is None after deserialization; construct directly
    (same approach as Phase 2 prototype make_modeller())."""
    from dials.algorithms.profile_model.gaussian_rs import GaussianRSProfileModeller
    from dials.array_family import flex as _flex

    GridMethod = GaussianRSProfileModeller.GridMethod
    FitMethod = GaussianRSProfileModeller.FitMethod
    grid_method = int(GridMethod.names["regular_grid"].real)
    fit_method = int(FitMethod.names["cell_cache_scatter"].real)
    profile = expr.profile
    if profile._scan_varying:
        sigma_b = _flex.mean(profile.sigma_b(deg=False))
        sigma_m = _flex.mean(profile.sigma_m(deg=False))
    else:
        sigma_b = profile.sigma_b(deg=False)
        sigma_m = profile.sigma_m(deg=False)
    return GaussianRSProfileModeller(
        expr.beam,
        expr.detector,
        expr.goniometer,
        expr.scan,
        sigma_b,
        sigma_m,
        profile.n_sigma() * 1.5,
        5,
        num_scan_points,
        0.02,
        grid_method,
        fit_method,
    )


class ChunkDriver:
    """Drives the single-pass integration loop.

    Each frame is read exactly once; each reflection processed exactly once.
    Slice finalization is triggered eagerly as the frame cursor advances.
    MVP: single experiment.
    """

    SLICE_WIDTH_FRAMES: int = 5  # matches production scan_step=5 default
    SLICES_PER_CHUNK: int = (
        5  # must cover max bbox half-width (~11 frames = 2.2 slices) + fanout
    )

    def __init__(self, experiments, master_reflections, params) -> None:
        assert len(experiments) == 1, "single-pass MVP: single experiment only"
        self.experiments = experiments
        self.master = master_reflections
        self.params = params

        scan = experiments[0].scan
        self._scan_lo, self._scan_hi = scan.get_array_range()
        sw = self.SLICE_WIDTH_FRAMES
        total_frames = self._scan_hi - self._scan_lo
        num_scan_points = (total_frames + sw - 1) // sw
        self._num_scan_points = num_scan_points

        self._sampler = _build_sampler(experiments, num_scan_points)
        self._profile_fitter = self._build_validated_modeller(
            experiments, num_scan_points
        )
        self._readiness_tracker = SliceReadinessTracker(num_scan_points, 3 * 3)
        self._pending_queue = PendingFitQueue(cells_per_slice=3 * 3)
        self._next_finalize = 0

    def _build_validated_modeller(self, experiments, num_scan_points: int):
        from dials.algorithms.integration.validation import (
            ValidatedMultiExpProfileModeller,
        )
        from dials.algorithms.profile_model.modeller import MultiExpProfileModeller

        pf = MultiExpProfileModeller()
        for expr in experiments:
            pf.add(_make_modeller(expr, num_scan_points))
        return ValidatedMultiExpProfileModeller([pf])

    def run(self):
        from dials.algorithms.integration.processor import ShoeboxProcessor
        from dials.array_family import flex
        from dials.model.data import make_image

        imageset = self.experiments[0].imageset
        scan_lo, scan_hi = self._scan_lo, self._scan_hi
        n_panels = len(imageset.get_detector())

        # Compute partiality — required by profile modeller's check1 filter
        self.master.compute_partiality(self.experiments)

        # Pre-create fitting output columns so pending-queue drain_by_slice
        # can write back via set_selected (Python KeyError if column absent).
        # The C++ set_selected_rows_index auto-creates for inline-fit, but
        # the Python path in PendingFitQueue needs them to exist already.
        n = len(self.master)
        for col in (
            "intensity.prf.value",
            "intensity.prf.variance",
            "profile.correlation",
        ):
            if col not in self.master:
                self.master[col] = flex.double(n, 0.0)

        # Tag each master row with its position so pending-queue writes can
        # target the correct master row after ShoeboxProcessor sub-tables
        # have been round-tripped through select_rows_index/set_selected_rows_index.
        self.master["_master_idx"] = flex.int(range(len(self.master)))

        # Single executor for the full scan
        executor = SinglePassExecutor(
            self.experiments,
            self._profile_fitter,
            self._sampler,
            self._readiness_tracker,
            self._pending_queue,
        )

        # Initialize and allocate shoeboxes for the full scan
        executor.initialize(scan_lo, scan_hi, self.master)
        self.master["shoebox"] = flex.shoebox(
            self.master["panel"],
            self.master["bbox"],
            allocate=False,
            flatten=False,
        )
        proc = ShoeboxProcessor(self.master, n_panels, scan_lo, scan_hi, False)

        # Process all frames sequentially — each frame read exactly once
        total = scan_hi - scan_lo
        log_interval = max(1, total // 10)
        for i in range(total):
            image = imageset.get_corrected_data(i)
            mask = imageset.get_mask(i)
            proc.next(make_image(image, mask), executor)
            del image, mask
            # Eagerly finalize slices whose contributions are complete
            self._maybe_finalize_slices(scan_lo + i + 1)
            if (i + 1) % log_interval == 0 or i + 1 == total:
                logger.info(" Processed %d/%d frames", i + 1, total)

        assert proc.finished(), "ShoeboxProcessor did not finish"
        del self.master["shoebox"]
        executor.finalize()

        # Finalize remaining slices (tail of scan) and flush pending queue
        self._finalize_remaining_slices()
        self._pending_queue.flush_all(self._profile_fitter.modellers[0], self.master)

        # Reconcile boundary reflections with baseline behavior.
        # fit_reciprocal_space() sets var=-1.0 as default before attempting fit;
        # boundary reflections throw dials::error in GridSampler::nearest(),
        # leaving var=-1.0 and setting failed_during_profile_fitting.
        # Our _handle_fit() skips these entirely, so replicate the sentinel here.
        prf_var = self.master["intensity.prf.variance"]
        prf_flag = self.master.get_flags(self.master.flags.integrated_prf)
        never_fitted = (prf_var == 0.0) & ~prf_flag
        prf_var.set_selected(never_fitted, -1.0)
        self.master.set_flags(
            never_fitted, self.master.flags.failed_during_profile_fitting
        )

        # Clean up the internal tracking column
        del self.master["_master_idx"]

        logger.info(
            " Single-pass complete: %d frames, %d reflections",
            scan_hi - scan_lo,
            len(self.master),
        )
        return self.master

    def _maybe_finalize_slices(self, current_frame):
        """Eagerly finalize slices whose modelling contributions are all complete.

        Slice k needs contributions from reflections with nearest() in slices
        {k-1, k, k+1}. Conservative trigger: finalize slice k when
        current_frame >= scan_lo + sw * (k + SLICES_PER_CHUNK).
        +3 accounts for: reflections in slice k+1 have xyzcal up to 5(k+2),
        bboxes extend ~5 frames max, completing by frame 5(k+3).
        """
        sw = self.SLICE_WIDTH_FRAMES
        max_safe_k = (current_frame - self._scan_lo) // sw - self.SLICES_PER_CHUNK
        while (
            self._next_finalize <= max_safe_k
            and self._next_finalize < self._num_scan_points
        ):
            self._finalize_slice(self._next_finalize)
            self._next_finalize += 1

    def _finalize_remaining_slices(self):
        """Finalize all not-yet-finalized slices (tail of scan)."""
        while self._next_finalize < self._num_scan_points:
            self._finalize_slice(self._next_finalize)
            self._next_finalize += 1

    def _finalize_slice(self, slice_idx: int) -> None:
        """Finalize all cells in the slice, then drain the pending queue.

        ValidatedMultiExpProfileModeller.modellers[0] is the MultiExpProfileModeller.
        MultiExpProfileModeller[0] is the GaussianRSProfileModeller for experiment 0.
        finalize_cell(j) normalizes data_[j] in-place; does NOT set the global
        finalized_ flag, so remaining cells stay writable.
        """
        cells_per_slice = self._readiness_tracker.cells_per_slice
        cell_lo = cells_per_slice * slice_idx
        cell_hi = cells_per_slice * (slice_idx + 1)

        # Access the underlying GaussianRSProfileModeller
        modeller = self._profile_fitter.modellers[0][0]

        for cell_idx in range(cell_lo, min(cell_hi, modeller.size())):
            modeller.finalize_cell(cell_idx)

        self._readiness_tracker.mark_slice_finalized(slice_idx)

        # Drain pending reflections for this slice — fit them now that cells
        # are finalized.  Pass the raw MultiExpProfileModeller so that
        # compute_fitted_intensity can call .fit() on it.
        self._pending_queue.drain_by_slice(
            slice_idx, self._profile_fitter.modellers[0], self.master
        )

        logger.debug(
            " Finalized slice %d (cells %d-%d)",
            slice_idx,
            cell_lo,
            cell_hi - 1,
        )


def run_single_pass_parallel(experiments, reflections, params, nproc):
    """Run single-pass integration in parallel with block splitting.

    Computes overlapping frame ranges for each worker, dispatches
    SinglePassBlockTask instances via multiprocessing, and assembles
    the combined reflections and profile modeller.

    Args:
        experiments: ExperimentList (single experiment)
        reflections: reflection_table with bbox, xyzcal.px, etc.
        params: Parameters object (from integrator.py)
        nproc: int, number of worker processes

    Returns:
        (combined_reflections, combined_modeller)
        where combined_modeller is a ValidatedMultiExpProfileModeller
        accumulated from all workers.
    """
    from dials.algorithms.integration.processor import execute_parallel_task
    from dials.util.log import rehandle_cached_records
    from dials.util.mp import multi_node_parallel_map

    # 1. Extract scan parameters
    scan = experiments[0].scan
    scan_lo, scan_hi = scan.get_array_range()
    total_frames = scan_hi - scan_lo

    # 2. Compute slice geometry (same as ChunkDriver)
    sw = ChunkDriver.SLICE_WIDTH_FRAMES
    num_scan_points = (total_frames + sw - 1) // sw

    # 3. Compute overlap
    # Model overlap: one full z-slice (GridSampler's +/-1 z-fanout)
    model_overlap = sw

    # Bbox overlap: 99th percentile of bbox z-widths
    bbox = reflections["bbox"]
    z0_col, z1_col = bbox.parts()[4], bbox.parts()[5]
    frame_widths = z1_col - z0_col
    sorted_widths = sorted(frame_widths)
    bbox_overlap = int(sorted_widths[int(0.99 * len(sorted_widths))])

    overlap = model_overlap + bbox_overlap

    # Clamp nproc if scan is too short; return None to signal caller
    # should fall back to the serial ChunkDriver path
    if nproc * 2 > total_frames:
        nproc = max(1, total_frames // 2)
    if nproc < 2:
        return None, None

    # 4. Compute core ranges (evenly divide total_frames)
    core_size = total_frames // nproc
    remainder = total_frames % nproc
    core_ranges = []
    start = scan_lo
    for i in range(nproc):
        size = core_size + (1 if i < remainder else 0)
        core_ranges.append((start, start + size))
        start += size

    # 5. Compute extended ranges (core +/- overlap, clipped to scan bounds)
    extended_ranges = []
    for core_lo, core_hi in core_ranges:
        ext_lo = max(scan_lo, core_lo - overlap)
        ext_hi = min(scan_hi, core_hi + overlap)
        extended_ranges.append((ext_lo, ext_hi))

    # 6. Compute partiality (needed before worker selection)
    reflections.compute_partiality(experiments)

    # 7. Select reflections per worker and create tasks
    tasks = []
    bbox = reflections["bbox"]
    z0_col, z1_col = bbox.parts()[4], bbox.parts()[5]
    for i in range(nproc):
        ext_lo, ext_hi = extended_ranges[i]
        core_lo, core_hi = core_ranges[i]
        # Strict containment: bbox fully within extended range
        sel = (z0_col >= ext_lo) & (z1_col <= ext_hi)
        worker_refls = reflections.select(sel)
        if len(worker_refls) == 0:
            continue
        task = SinglePassBlockTask(
            index=i,
            experiments=experiments,
            reflections=worker_refls,
            extended_range=(ext_lo, ext_hi),
            core_range=(core_lo, core_hi),
            num_scan_points=num_scan_points,
        )
        tasks.append(task)

    logger.info(
        " Parallel single-pass: %d workers, %d frames, overlap=%d",
        len(tasks),
        total_frames,
        overlap,
    )

    # 8. Dispatch via multi_node_parallel_map
    results = []
    n_tasks = len(tasks)

    def process_output(result_pair):
        rehandle_cached_records(result_pair[1])
        results.append(result_pair[0])
        logger.info(" Worker %d/%d complete", len(results), n_tasks)

    multi_node_parallel_map(
        func=execute_parallel_task,
        iterable=tasks,
        njobs=1,
        nproc=nproc,
        callback=process_output,
        cluster_method="multiprocessing",
        preserve_order=True,
    )

    # 9. Assemble results
    if not results:
        return None, None
    results.sort(key=lambda r: r.index)

    # Combine reflections
    combined = results[0].reflections
    for r in results[1:]:
        combined.extend(r.reflections)

    # Accumulate modellers for report
    combined_modeller = None
    for r in results:
        if r.data is None:
            continue
        if combined_modeller is None:
            combined_modeller = r.data
        else:
            combined_modeller.accumulate(r.data)

    logger.info(
        " Parallel single-pass complete: %d reflections from %d workers",
        len(combined),
        len(results),
    )

    return combined, combined_modeller

"""Phase 3 MVP: single-pass chunked modeller-integrator skeleton."""

from __future__ import annotations

import logging

from dials.algorithms.shoebox import MaskCode
from dials_algorithms_integration_integrator_ext import Executor

logger = logging.getLogger(__name__)


class SliceReadinessTracker:
    """Tracks finalized profile-model z-slices.
    Task 3 stub: is_cell_ready() always returns False. Task 4 replaces."""

    def __init__(self, num_slices: int, cells_per_slice: int) -> None:
        self.num_slices = num_slices
        self.cells_per_slice = cells_per_slice
        self.finalized_slices: set[int] = set()

    def is_cell_ready(self, cell_idx: int) -> bool:
        return (cell_idx // self.cells_per_slice) in self.finalized_slices

    def mark_slice_finalized(self, slice_idx: int) -> None:
        self.finalized_slices.add(slice_idx)


class PendingFitQueue:
    """Holds reflections pending profile-model cell finalization.
    Task 3 stub: enqueue/drain are no-ops. Task 4 replaces."""

    def __init__(self) -> None:
        self._by_cell: dict[int, list] = {}  # cell_idx -> list of row entries

    def enqueue(self, reflections, nearest_cells) -> None:
        pass  # Task 4

    def drain_by_slice(self, slice_idx, profile_fitter, master_table) -> None:
        pass  # Task 4

    def flush_all(self, profile_fitter, master_table) -> None:
        pass  # Task 4

    def is_empty(self) -> bool:
        return len(self._by_cell) == 0


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

        # Step 11: fitting (Task 3 stub: all reflections, inline, unconditional)
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
        """Task 3 stub: fit all reflections inline unconditionally.
        Task 4 splits into inline-fit vs pending-queue based on cell readiness.
        Uses modellers[0] (the raw MultiExpProfileModeller) because
        compute_fitted_intensity requires an object with a .fit() method;
        ValidatedMultiExpProfileModeller delegates modelling but not fitting."""
        reflections.compute_fitted_intensity(self.profile_fitter.modellers[0])

    def finalize(self) -> None:
        pass

    def data(self):
        return self.profile_fitter


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
    SLICES_PER_CHUNK: int = 3

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
        self._pending_queue = PendingFitQueue()
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
        for i in range(scan_hi - scan_lo):
            image = imageset.get_corrected_data(i)
            mask = imageset.get_mask(i)
            proc.next(make_image(image, mask), executor)
            del image, mask
            # Eagerly finalize slices whose contributions are complete
            self._maybe_finalize_slices(scan_lo + i + 1)

        assert proc.finished(), "ShoeboxProcessor did not finish"
        del self.master["shoebox"]
        executor.finalize()

        # Finalize remaining slices (tail of scan) and flush pending queue
        self._finalize_remaining_slices()
        self._pending_queue.flush_all(self._profile_fitter, self.master)

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
        """Task 3 stub: mark slice finalized in tracker.
        Task 4 adds: finalize_cell() per cell, then drain pending_queue."""
        self._readiness_tracker.mark_slice_finalized(slice_idx)

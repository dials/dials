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


def _process_chunk(experiments, chunk_reflections, imageset_slice, executor):
    """Run ShoeboxProcessor over imageset_slice, dispatching to executor per frame.
    Lifted from Phase 2 prototype process_chunk() (lines 33-61)."""
    from dials.algorithms.integration.processor import ShoeboxProcessor
    from dials.array_family import flex
    from dials.model.data import make_image

    frame0, frame1 = imageset_slice.get_array_range()
    executor.initialize(frame0, frame1, chunk_reflections)
    chunk_reflections["shoebox"] = flex.shoebox(
        chunk_reflections["panel"],
        chunk_reflections["bbox"],
        allocate=False,
        flatten=False,
    )
    proc = ShoeboxProcessor(
        chunk_reflections, len(imageset_slice.get_detector()), frame0, frame1, False
    )
    for i in range(len(imageset_slice)):
        image = imageset_slice.get_corrected_data(i)
        mask = imageset_slice.get_mask(i)
        proc.next(make_image(image, mask), executor)
        del image, mask
    assert proc.finished(), "ShoeboxProcessor did not finish"
    del chunk_reflections["shoebox"]
    executor.finalize()
    return chunk_reflections, executor.data()


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
    """Drives the single-pass chunked integration loop.

    Chunk N covers nominal frames [5N, 5N+15) = z-slices {N, N+1, N+2}.
    After chunk N completes, z-slice N is finalized (Task 4 wires the real logic).
    Chunks advance 1 slice per iteration (sliding window). MVP: single experiment.
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
        """Run all chunks; return updated master reflection table."""
        from dials.array_family import flex

        scan_lo, sw, num_slices = (
            self._scan_lo,
            self.SLICE_WIDTH_FRAMES,
            self._num_scan_points,
        )
        imageset = self.experiments[0].imageset
        allowed_lo, _hi = imageset.get_array_range()
        total_dispatched = chunks_processed = 0

        for chunk_idx in range(num_slices):
            nominal_lo = scan_lo + sw * chunk_idx
            frame_lo = max(nominal_lo, self._scan_lo)
            frame_hi = min(nominal_lo + sw * self.SLICES_PER_CHUNK, self._scan_hi)
            if frame_lo >= frame_hi:
                continue

            bboxes = self.master["bbox"]
            bbox_lo = flex.int([bb[4] for bb in bboxes])
            bbox_hi = flex.int([bb[5] for bb in bboxes])
            chunk_reflections = self.master.select(
                (bbox_lo < frame_hi) & (bbox_hi > frame_lo)
            )

            if len(chunk_reflections) == 0:
                continue

            # Use actual bbox extents (not nominal range) as the imageset slice bounds.
            # ShoeboxProcessor asserts bbox[4] >= frame0, so frame0 must be the
            # minimum bbox start of all selected reflections.
            chunk_bboxes = chunk_reflections["bbox"]
            actual_lo = min(bb[4] for bb in chunk_bboxes)
            actual_hi = max(bb[5] for bb in chunk_bboxes)

            logger.info(
                " Chunk %d: frames [%d, %d), %d reflections",
                chunk_idx,
                actual_lo,
                actual_hi,
                len(chunk_reflections),
            )
            total_dispatched += len(chunk_reflections)
            imageset_slice = imageset[actual_lo - allowed_lo : actual_hi - allowed_lo]
            executor = SinglePassExecutor(
                self.experiments,
                self._profile_fitter,
                self._sampler,
                self._readiness_tracker,
                self._pending_queue,
            )
            _process_chunk(
                self.experiments, chunk_reflections, imageset_slice, executor
            )
            chunks_processed += 1
            self._finalize_slice(chunk_idx)

        # Post-loop: finalize trailing slices; flush pending queue (stubs in Task 3)
        for slice_idx in range(
            max(0, num_slices - self.SLICES_PER_CHUNK + 1), num_slices
        ):
            self._finalize_slice(slice_idx)
        self._pending_queue.flush_all(self._profile_fitter, self.master)
        logger.info(
            " Single-pass complete: %d chunks, %d reflection-frame dispatches",
            chunks_processed,
            total_dispatched,
        )
        return self.master

    def _finalize_slice(self, slice_idx: int) -> None:
        """Task 3 stub: mark slice finalized in tracker.
        Task 4 adds: finalize_cell() per cell, then drain pending_queue."""
        self._readiness_tracker.mark_slice_finalized(slice_idx)

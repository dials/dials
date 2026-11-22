from __future__ import annotations

import math
import os
from unittest import mock

import numpy as np
import psutil
import pytest

from dxtbx.model.experiment_list import ExperimentListFactory

import dials.algorithms.integration.processor
from dials.algorithms.profile_model.gaussian_rs import Model
from dials.array_family import flex
from dials.model.data import make_image
from dials_algorithms_integration_integrator_ext import (
    Executor,
    JobList,
    ShoeboxProcessor,
    max_memory_needed,
)


def test_shoebox_memory_is_a_reasonable_guesstimate(dials_data):
    path = dials_data("centroid_test_data") / "experiments.json"

    exlist = ExperimentListFactory.from_json_file(path)[0]
    exlist.profile = Model(
        None,
        n_sigma=3,
        sigma_b=0.024 * math.pi / 180.0,
        sigma_m=0.044 * math.pi / 180.0,
    )

    rlist = flex.reflection_table.from_predictions(exlist)
    rlist["id"] = flex.int(len(rlist), 0)
    rlist["bbox"] = flex.int6(rlist.size(), (0, 1, 0, 1, 0, 1))

    jobs = JobList()
    jobs.add((0, 1), (0, 9), 9, 0)
    # Every single-pixel shoebox ends on frame 0, so all 1996 of them hold
    # data, mask and background (9 bytes per pixel) at the same time, on top
    # of the Shoebox object and one frame lookup index (8 bytes) each
    expected = 1996 * (9 + 8 + _shoebox_object_bytes())
    for flatten in (True, False):
        assumed_memory_usage = list(jobs.shoebox_memory(rlist, flatten))
        assert assumed_memory_usage == [expected]

        max_mem = max_memory_needed(rlist, 0, 9, flatten)
        assert max_mem == expected


def _reflections(bboxes):
    rlist = flex.reflection_table()
    rlist["id"] = flex.int(len(bboxes), 0)
    rlist["flags"] = flex.size_t(len(bboxes), 0)
    rlist["bbox"] = flex.int6(bboxes)
    return rlist


def _shoebox_object_bytes():
    """The fixed cost of an unallocated Shoebox object, which is derived
    from sizeof() so depends on the platform"""
    # A single pixel on a single frame costs 9 bytes of pixel data, plus one
    # frame lookup index
    return max_memory_needed(_reflections([(0, 1, 0, 1, 0, 1)]), 0, 1, False) - 17


def test_shoebox_object_bytes_is_plausible():
    # Measured as 448 bytes per unallocated shoebox on 64-bit Linux
    assert 300 < _shoebox_object_bytes() < 600


@pytest.mark.parametrize(
    "flatten,expected_pixel_memory",
    [
        # Peak is on the last frame of A: data and mask (5 bytes per pixel) for
        # all 4 frames of A, plus background (4 bytes per pixel). B has been freed.
        (False, 5 * 400 + 4 * 400),
        # Flattened shoeboxes have a single frame of 100 pixels for A
        (True, 5 * 100 + 4 * 100),
    ],
)
def test_shoebox_memory_counts_background_only_on_final_frame(
    flatten, expected_pixel_memory
):
    rlist = _reflections(
        [
            (0, 10, 0, 10, 0, 4),  # A: 10x10 pixels on frames 0-3
            (0, 5, 0, 5, 1, 2),  # B: 5x5 pixels on frame 1
        ]
    )
    # Both Shoebox objects and the frame lookup (4 + 1 entries) exist for the
    # whole job
    expected = expected_pixel_memory + 2 * _shoebox_object_bytes() + 5 * 8

    jobs = JobList()
    jobs.add((0, 1), (0, 4), 4, 0)
    assert list(jobs.shoebox_memory(rlist, flatten)) == [expected]
    assert max_memory_needed(rlist, 0, 4, flatten) == expected


def test_shoebox_memory_includes_bytes_per_reflection():
    rlist = _reflections([(0, 10, 0, 10, 0, 4), (0, 5, 0, 5, 1, 2)])
    jobs = JobList()
    jobs.add((0, 1), (0, 4), 4, 0)
    for flatten in (True, False):
        base = max_memory_needed(rlist, 0, 4, flatten)
        assert max_memory_needed(rlist, 0, 4, flatten, 100) == base + 2 * 100
        assert list(jobs.shoebox_memory(rlist, flatten, 100)) == [base + 2 * 100]


def test_shoebox_memory_with_compression():
    rlist = _reflections([(0, 10, 0, 10, 0, 4), (0, 5, 0, 5, 1, 2)])
    jobs = JobList()
    jobs.add((0, 1), (0, 4), 4, 0)
    # The Shoebox objects and frame lookup as before, plus the per-reflection
    # record of the compressed store
    fixed = 2 * (_shoebox_object_bytes() + 12) + 5 * 8
    # Peak is on the last frame of A: frames 0-2 of A compressed, at an assumed
    # 2 bytes per pixel plus one 256-byte chunk, and all of A decoded with its
    # background (9 bytes per pixel). B is on one frame, so never compressed.
    expected = fixed + 2 * 300 + 256 + 9 * 400
    assert max_memory_needed(rlist, 0, 4, False, 0, True) == expected
    assert list(jobs.shoebox_memory(rlist, False, 0, True)) == [expected]
    # Flat shoeboxes are a single slice, so never compressed
    assert max_memory_needed(rlist, 0, 4, True, 0, True) == max_memory_needed(
        rlist, 0, 4, True
    )


class _RecordingExecutor(Executor):
    """Keep a copy of each shoebox as the executor sees it"""

    def __init__(self):
        super().__init__()
        self.seen = {}
        self.calls = []

    def process(self, frame, reflections):
        self.calls.append((frame, len(reflections)))
        for i, sbox in zip(reflections["index"], reflections["shoebox"]):
            assert sbox.is_consistent()
            self.seen[i] = (
                frame,
                sbox.data.as_numpy_array().copy(),
                sbox.mask.as_numpy_array().copy(),
                sbox.background.as_numpy_array().copy(),
            )


def _extract(images, masks, bboxes, compress, save=False):
    rlist = flex.reflection_table()
    rlist["bbox"] = flex.int6(bboxes)
    rlist["panel"] = flex.size_t(len(bboxes), 0)
    rlist["index"] = flex.size_t_range(len(bboxes))
    rlist["shoebox"] = flex.shoebox(
        rlist["panel"], rlist["bbox"], allocate=False, flatten=False
    )
    processor = ShoeboxProcessor(rlist, 1, 0, len(images), save, compress)
    assert processor.compressed() == compress
    executor = _RecordingExecutor()
    for image, mask in zip(images, masks):
        processor.next(make_image(image, mask), executor)
    assert processor.finished()
    return executor, rlist, processor


@pytest.mark.parametrize("save", [False, True])
def test_shoebox_processor_compression_is_lossless(save):
    rng = np.random.default_rng(42)
    ny, nx, nz = 40, 50, 6
    # Poisson counts on some frames and arbitrary floats on others, so that
    # both compressible and incompressible slices are stored
    images, masks = [], []
    for z in range(nz):
        if z % 2:
            data = rng.normal(0, 1e6, (ny, nx))
        else:
            data = rng.poisson(3, (ny, nx)).astype(float)
        images.append((flex.double(data),))
        masks.append((flex.bool(rng.random((ny, nx)) > 0.1),))
    bboxes = [
        (5, 15, 5, 20, 0, 6),  # whole scan
        (10, 12, 30, 31, 2, 3),  # one frame
        (-4, 6, -3, 8, 1, 5),  # off the image at the top left
        (45, 58, 35, 44, 0, 4),  # off the image at the bottom right
        (20, 30, 10, 20, 3, 6),  # ends on the last frame
        (0, 50, 0, 40, 1, 3),  # the whole image
    ]
    executor_c, rlist_c, processor = _extract(images, masks, bboxes, True, save)
    executor_u, rlist_u, _ = _extract(images, masks, bboxes, False, save)
    seen_c, seen_u = executor_c.seen, executor_u.seen
    assert sorted(seen_c) == sorted(seen_u) == list(range(len(bboxes)))
    for i in seen_u:
        assert seen_c[i][0] == seen_u[i][0]
        for a, b in zip(seen_c[i][1:], seen_u[i][1:]):
            assert a.dtype == b.dtype and np.array_equal(a, b)
    for sc, su in zip(rlist_c["shoebox"], rlist_u["shoebox"]):
        assert sc.is_data_allocated() == su.is_data_allocated() == save
        if save:
            assert sc.data.all_eq(su.data) and sc.mask.all_eq(su.mask)
    assert processor.max_compressed_bytes() > 0


def test_shoebox_processor_compresses_counts():
    rng = np.random.default_rng(0)
    ny, nx, nz = 100, 100, 5
    images = [(flex.double(rng.poisson(1, (ny, nx)).astype(float)),) for _ in range(nz)]
    masks = [(flex.bool(flex.grid(ny, nx), True),) for _ in range(nz)]
    bboxes = [
        (x, x + 20, y, y + 20, 0, nz)
        for x in range(0, 80, 20)
        for y in range(0, 80, 20)
    ]
    _, _, processor = _extract(images, masks, bboxes, True)
    # 16 shoeboxes of 4 stored 20x20 slices: raw data and mask take 5 bytes
    # per pixel, but low counts should compress to under 1 byte per pixel
    raw = 16 * 4 * 400 * 5
    assert processor.max_compressed_bytes() < raw / 5


def test_shoebox_processor_batches_frames_where_many_shoeboxes_finish():
    # Five 1000x1000x2 shoeboxes finish on the last frame: 18 MB each decoded,
    # 90 MB in all, far more than on a typical frame and more than the 64 MB
    # lower limit for one batch. So they are decoded and processed in two
    # batches, with the same results.
    rng = np.random.default_rng(1)
    ny, nx, nz = 1000, 1000, 3
    images = [(flex.double(rng.poisson(2, (ny, nx)).astype(float)),) for _ in range(nz)]
    masks = [(flex.bool(rng.random((ny, nx)) > 0.01),) for _ in range(nz)]
    bboxes = [(0, nx, 0, ny, 1, 3)] * 5 + [(0, 5, 0, 5, z, z + 1) for z in range(3)]
    executor_c, _, _ = _extract(images, masks, bboxes, True)
    executor_u, _, _ = _extract(images, masks, bboxes, False)
    assert executor_u.calls == [(0, 1), (1, 1), (2, 6)]
    assert executor_c.calls == [(0, 1), (1, 1), (2, 3), (2, 3)]
    for i in executor_u.seen:
        for a, b in zip(executor_c.seen[i][1:], executor_u.seen[i][1:]):
            assert np.array_equal(a, b)


def test_shoebox_processor_never_compresses_flat_shoeboxes():
    rlist = flex.reflection_table()
    rlist["bbox"] = flex.int6([(0, 2, 0, 2, 0, 3)])
    rlist["panel"] = flex.size_t(1, 0)
    rlist["shoebox"] = flex.shoebox(
        rlist["panel"], rlist["bbox"], allocate=False, flatten=True
    )
    assert not ShoeboxProcessor(rlist, 1, 0, 3, False, True).compressed()


@mock.patch("dials.algorithms.integration.processor.flex.max")
@mock.patch("dials.algorithms.integration.processor.MEMORY_LIMIT", 1_000_000_000)
def test_runtime_error_raised_when_not_enough_memory(mock_flex_max):
    # We are taking into account current memory usage, which will typically be ~200MB.
    # So set the shoebox memory limit to be 500MB, which should fail the memory test.
    mock_flex_max.return_value = 500_000_000

    phil_mock = mock.Mock()
    phil_mock.mp.method = "multiprocessing"
    phil_mock.mp.nproc = 4
    phil_mock.block.max_memory_usage = 0.5

    reflections = {"bbox": flex.int6(1000, (0, 1, 0, 1, 0, 1))}
    # No detectors, so no memory is required for images
    experiments = mock.Mock()
    experiments.detectors.return_value = []
    manager = dials.algorithms.integration.processor._Manager(
        experiments, reflections, phil_mock
    )
    manager.jobs = mock.Mock(autospec=JobList)

    with pytest.raises(MemoryError) as exc_info:
        manager.compute_processors()
    assert "Not enough memory to run integration jobs." in exc_info.value.args[0]
    mock_flex_max.assert_called_once_with(manager.jobs.shoebox_memory.return_value)

    # Reduce memory usage of the shoeboxes by the current memory usage
    # (with an extra 1MB for safety), should then pass.
    current_memory_usage = psutil.Process(os.getpid()).memory_info().rss
    mock_flex_max.return_value = 499_000_000 - current_memory_usage
    manager.compute_processors()
    mock_flex_max.assert_called_with(manager.jobs.shoebox_memory.return_value)

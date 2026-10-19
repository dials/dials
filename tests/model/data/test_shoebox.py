from __future__ import annotations

import math
import pickle
import random

import pytest

from scitbx import matrix

from dials.array_family import flex
from dials.model.data import FrameSlicedShoebox, Shoebox


def random_shoeboxes(num, mask=False):
    for i in range(num):
        x0 = random.randint(0, 100)
        y0 = random.randint(0, 100)
        z0 = random.randint(0, 100)
        x1 = random.randint(x0 + 5, x0 + 20)
        y1 = random.randint(y0 + 5, y0 + 20)
        z1 = random.randint(z0 + 5, z0 + 20)
        bbox = (x0, x1, y0, y1, z0, z1)
        xc0 = (x1 + x0) / 2.0
        yc0 = (y1 + y0) / 2.0
        zc0 = (z1 + z0) / 2.0
        xc = random.uniform(xc0 - 1, xc0 + 1)
        yc = random.uniform(yc0 - 1, yc0 + 1)
        zc = random.uniform(zc0 - 1, zc0 + 1)
        centre = (xc, yc, zc)
        intensity = random.randint(10, 10000)
        shoebox = generate_shoebox(bbox, centre, intensity, mask=mask)
        yield (shoebox, (centre, intensity))


def generate_shoebox(bbox, centre, intensity, mask=False):
    from dials.algorithms.shoebox import MaskCode

    shoebox = Shoebox()
    shoebox.bbox = bbox
    shoebox.allocate_data()
    shoebox.allocate_background()
    for i in range(len(shoebox.mask)):
        shoebox.mask[i] = MaskCode.Valid | MaskCode.Foreground
    shoebox.data = gaussian(
        shoebox.size(),
        1.0,
        [c - o for c, o in zip(centre[::-1], shoebox.offset())],
        [s / 8.0 for s in shoebox.size()],
    )
    if mask:
        shoebox.mask = create_mask(
            shoebox.size(),
            [c - o for c, o in zip(centre[::-1], shoebox.offset())],
            MaskCode.Valid | MaskCode.Foreground,
        )
    tot = 0
    mask_code = MaskCode.Valid | MaskCode.Foreground
    for i in range(len(shoebox.data)):
        if shoebox.mask[i] & mask_code == mask_code:
            tot += shoebox.data[i]
    if tot > 0:
        shoebox.data *= intensity / tot
    return shoebox


def create_mask(size, x0, value):
    from scitbx.array_family import flex

    mask = flex.uint8(flex.grid(size), 0)
    rad = min(s - c for s, c in zip(size, x0))
    for k in range(size[0]):
        for j in range(size[1]):
            for i in range(size[2]):
                d = math.sqrt((j - x0[1]) ** 2 + (i - x0[2]) ** 2)
                if d < rad:
                    mask[k, j, i] = value
    return mask


def evaluate_gaussian(x, a, x0, sx):
    assert len(x) == len(x0)
    assert len(x) == len(sx)

    g = 0.0
    for xi, x0i, sxi in zip(x, x0, sx):
        g += (xi - x0i) ** 2 / (2.0 * sxi**2)

    return a * math.exp(-g)


def gaussian(size, a, x0, sx):
    from dials.array_family import flex

    result = flex.real(flex.grid(size))

    index = [0] * len(size)
    while True:
        result[index] = evaluate_gaussian(index, a, x0, sx)
        for j in range(len(size)):
            index[j] += 1
            if index[j] < size[j]:
                break
            index[j] = 0
            if j == len(size) - 1:
                return result


def test_allocate():
    for i in range(10):
        x0 = random.randint(0, 1000)
        y0 = random.randint(0, 1000)
        z0 = random.randint(0, 1000)
        x1 = random.randint(1, 10) + x0
        y1 = random.randint(1, 10) + y0
        z1 = random.randint(1, 10) + z0

        shoebox = Shoebox((x0, x1, y0, y1, z0, z1))
        shoebox.allocate_data()
        shoebox.allocate_background()
        assert shoebox.data.all() == (z1 - z0, y1 - y0, x1 - x0)
        assert shoebox.mask.all() == (z1 - z0, y1 - y0, x1 - x0)
        shoebox.deallocate()
        assert shoebox.data.all() == (0, 0, 0)
        assert shoebox.mask.all() == (0, 0, 0)


def test_offset():
    for i in range(10):
        x0 = random.randint(0, 1000)
        y0 = random.randint(0, 1000)
        z0 = random.randint(0, 1000)
        x1 = random.randint(1, 10) + x0
        y1 = random.randint(1, 10) + y0
        z1 = random.randint(1, 10) + z0

        shoebox = Shoebox((x0, x1, y0, y1, z0, z1))
        assert shoebox.xoffset() == x0
        assert shoebox.yoffset() == y0
        assert shoebox.zoffset() == z0
        assert shoebox.offset() == (z0, y0, x0)


def test_size():
    for i in range(10):
        x0 = random.randint(0, 1000)
        y0 = random.randint(0, 1000)
        z0 = random.randint(0, 1000)
        x1 = random.randint(1, 10) + x0
        y1 = random.randint(1, 10) + y0
        z1 = random.randint(1, 10) + z0

        shoebox = Shoebox((x0, x1, y0, y1, z0, z1))
        assert shoebox.xsize() == x1 - x0
        assert shoebox.ysize() == y1 - y0
        assert shoebox.zsize() == z1 - z0
        assert shoebox.size() == (z1 - z0, y1 - y0, x1 - x0)


def test_consistent():
    from dials.array_family import flex

    for i in range(1000):
        x0 = random.randint(0, 1000)
        y0 = random.randint(0, 1000)
        z0 = random.randint(0, 1000)
        x1 = random.randint(1, 10) + x0
        y1 = random.randint(1, 10) + y0
        z1 = random.randint(1, 10) + z0
        try:
            shoebox = Shoebox((x0, x1, y0, y1, z0, z1))
            assert not shoebox.is_consistent()
            shoebox.allocate_data()
            shoebox.allocate_background()
            assert shoebox.is_consistent()
            shoebox.data = flex.real(flex.grid(20, 20, 20))
            assert not shoebox.is_consistent()
            shoebox.deallocate()
            assert not shoebox.is_consistent()
        except Exception:
            print(x0, y0, z0, x1, y1, z1)
            raise


def test_is_bbox_within_image_volume():
    isize = (1000, 1000)
    srange = (0, 100)

    shoebox = Shoebox((10, 20, 10, 20, 10, 20))
    assert shoebox.is_bbox_within_image_volume(isize, srange)
    shoebox = Shoebox((-10, 20, 10, 20, 10, 20))
    assert not shoebox.is_bbox_within_image_volume(isize, srange)
    shoebox = Shoebox((10, 20, -10, 20, 10, 20))
    assert not shoebox.is_bbox_within_image_volume(isize, srange)
    shoebox = Shoebox((10, 20, 10, 20, -10, 20))
    assert not shoebox.is_bbox_within_image_volume(isize, srange)
    shoebox = Shoebox((10, 1020, 10, 20, 10, 20))
    assert not shoebox.is_bbox_within_image_volume(isize, srange)
    shoebox = Shoebox((10, 20, 10, 1020, 10, 20))
    assert not shoebox.is_bbox_within_image_volume(isize, srange)
    shoebox = Shoebox((10, 20, 10, 20, 10, 1020))
    assert not shoebox.is_bbox_within_image_volume(isize, srange)


def test_does_bbox_contain_bad_pixels():
    from scitbx.array_family import flex

    mask = flex.bool(flex.grid(100, 100), True)
    for j in range(100):
        for i in range(40, 60):
            mask[j, i] = False
            mask[i, j] = False

    for i in range(1000):
        x0 = random.randint(0, 90)
        y0 = random.randint(0, 90)
        z0 = random.randint(0, 90)
        x1 = random.randint(1, 10) + x0
        y1 = random.randint(1, 10) + y0
        z1 = random.randint(1, 10) + z0

        shoebox = Shoebox((x0, x1, y0, y1, z0, z1))

        res1 = shoebox.does_bbox_contain_bad_pixels(mask)
        res2 = False
        if x0 >= 40 and x0 < 60:
            res2 = True
        if x1 > 40 and x1 <= 60:
            res2 = True
        if y0 >= 40 and y0 < 60:
            res2 = True
        if y1 > 40 and y1 <= 60:
            res2 = True

        assert res1 == res2


def test_count_mask_values():
    for i in range(10):
        x0 = random.randint(0, 90)
        y0 = random.randint(0, 90)
        z0 = random.randint(0, 90)
        x1 = random.randint(1, 10) + x0
        y1 = random.randint(1, 10) + y0
        z1 = random.randint(1, 10) + z0

        shoebox = Shoebox((x0, x1, y0, y1, z0, z1))
        shoebox.allocate_data()
        shoebox.allocate_background()
        maxnum = len(shoebox.mask)
        num = random.randint(1, maxnum)
        indices = random.sample(list(range(maxnum)), num)
        value = 1 << 2
        for i in indices:
            shoebox.mask[i] = value

        assert shoebox.count_mask_values(value) == num


def test_centroid_all():
    for shoebox, (XC, I) in random_shoeboxes(10):
        centroid = shoebox.centroid_all()
        assert shoebox.is_consistent()
        assert abs(matrix.col(centroid.px.position) - matrix.col(XC)) < 1.0


def test_centroid_masked():
    for shoebox, (XC, I) in random_shoeboxes(10):
        centroid = shoebox.centroid_masked(1 << 0)
        assert shoebox.is_consistent()
        assert abs(matrix.col(centroid.px.position) - matrix.col(XC)) < 1.0


def test_summed_intensity():
    for shoebox, (XC, I) in random_shoeboxes(10):
        intensity = shoebox.summed_intensity()
        assert shoebox.is_consistent()
        assert abs(intensity.observed.value - I) < 1e-1


def test_flatten():
    from dials.algorithms.shoebox import MaskCode
    from dials.array_family import flex

    for shoebox, (XC, I) in random_shoeboxes(10, mask=True):
        assert not shoebox.flat
        zs = shoebox.zsize()
        ys = shoebox.ysize()
        xs = shoebox.xsize()
        expected_data = flex.real(flex.grid(1, ys, xs), 0)
        expected_mask = flex.uint8(flex.grid(1, ys, xs), 0)
        for k in range(zs):
            for j in range(ys):
                for i in range(xs):
                    expected_data[0, j, i] += shoebox.data[k, j, i]
                    expected_mask[0, j, i] |= shoebox.mask[k, j, i]
                    if not (expected_mask[0, j, i] & MaskCode.Valid) or not (
                        shoebox.mask[k, j, i] & MaskCode.Valid
                    ):
                        expected_mask[0, j, i] &= ~MaskCode.Valid
        shoebox.flatten()
        diff = expected_data.as_double() - shoebox.data.as_double()
        max_diff = flex.max(flex.abs(diff))
        assert max_diff < 1e-7
        assert expected_mask.all_eq(shoebox.mask)
        assert shoebox.flat
        assert shoebox.is_consistent()


def test_all_foreground_valid():
    from .all_foreground_valid_data import data

    shoeboxes = pickle.loads(bytes(data, encoding="latin-1"), encoding="bytes")
    for i, shoebox in enumerate(shoeboxes):
        if i < 4:
            assert not shoebox.all_foreground_valid()
        else:
            assert shoebox.all_foreground_valid()


def frame_sliced_shoebox_test_data():
    """A 4 x 3 x 3 (x, y, z) shoebox, with the first frame at image 11"""
    from dials.algorithms.shoebox import MaskCode

    shoebox = Shoebox(0, (0, 4, 0, 3, 10, 13))
    shoebox.allocate_data()
    shoebox.allocate_background()

    # Every pixel on frame k has a value of k + 1 and a background of 0.5
    for k in range(3):
        for j in range(3):
            for i in range(4):
                shoebox.data[k, j, i] = k + 1
                shoebox.background[k, j, i] = 0.5
                shoebox.mask[k, j, i] = MaskCode.Valid | MaskCode.Background

    # Two foreground pixels on frame 0 and three on frame 1
    for k, j, i in [(0, 0, 0), (0, 0, 1), (1, 1, 0), (1, 1, 1), (1, 1, 2)]:
        shoebox.mask[k, j, i] = MaskCode.Valid | MaskCode.Foreground

    # A single foreground pixel on frame 2, which is not valid
    shoebox.mask[2, 2, 3] = MaskCode.Foreground

    return shoebox


def frame_sliced_shoebox_test_models(
    first_frame=11, nframes=3, phi_cal=1.15, sigma_phi=0.1
):
    """Everything the FrameSlicedShoebox constructor wants besides the shoebox
    and its Miller index, as keyword arguments. That is the rocking curve of the
    reflection, then the experimental models over the whole of a scan of nframes
    images starting at first_frame, alongside the image number that the first of
    those refers to.

    The scan has 0.1 radian oscillations beginning at zero, the crystal is
    stationary with a 10 Angstrom cubic cell aligned with the laboratory axes,
    and the beam is stationary along -z with a wavelength of 1 Angstrom"""
    return {
        "phi_cal": phi_cal,
        "sigma_phi": sigma_phi,
        "phi": flex.double(
            [0.1 * (f - 0.5) for f in range(first_frame, first_frame + nframes)]
        ),
        "phi_scan_points": flex.double(
            [0.1 * (f - 1) for f in range(first_frame, first_frame + nframes + 1)]
        ),
        "UB": flex.mat3_double([(0.1, 0, 0, 0, 0.1, 0, 0, 0, 0.1)] * nframes),
        "s0": flex.vec3_double([(0, 0, -1)] * nframes),
        "first_frame": first_frame,
    }


def gaussian_cdf(x):
    """The cumulative distribution function of the standard normal distribution"""
    return 0.5 * (1.0 + math.erf(x / math.sqrt(2.0)))


def test_frame_sliced_shoebox():
    sliced = FrameSlicedShoebox(
        frame_sliced_shoebox_test_data(),
        (1, 0, 0),
        **frame_sliced_shoebox_test_models(),
    )

    assert sliced.size() == 3
    assert len(sliced) == 3
    # The shoebox z range is 10 -> 13, which is images 11 -> 13 (one-based)
    assert list(sliced.frames) == [11, 12, 13]
    assert list(sliced.phi) == pytest.approx([1.05, 1.15, 1.25])
    # A stationary crystal and beam give the same excitation error on each frame
    assert list(sliced.excitation_error) == pytest.approx([math.sqrt(0.99) - 1.0] * 3)
    # The rocking curve is centred on the middle frame and is one frame wide, so
    # the middle frame takes the largest share of the reflection
    assert list(sliced.partiality) == pytest.approx(
        [
            gaussian_cdf(-0.5) - gaussian_cdf(-1.5),
            gaussian_cdf(0.5) - gaussian_cdf(-0.5),
            gaussian_cdf(1.5) - gaussian_cdf(0.5),
        ]
    )
    # No pixel is marked as used for the background, so the m/n term of the
    # variance is zero and the variance is just the sum of the raw foreground
    # pixel values
    assert list(sliced.summation_intensity) == [1.0, 4.5, 0.0]
    assert list(sliced.summation_intensity_variance) == [2.0, 6.0, 0.0]
    # The foreground pixel on frame 2 could not be summed, as it is not valid
    assert list(sliced.summation_intensity_valid) == [True, True, False]


def test_frame_sliced_shoebox_summation_intensity_counts_background_pixels():
    from dials.algorithms.shoebox import MaskCode

    # One frame, with two foreground pixels and four background pixels, of which
    # only two were used to calculate the background
    shoebox = Shoebox(0, (0, 6, 0, 1, 0, 1))
    shoebox.allocate_data()
    shoebox.allocate_background()
    for i in range(6):
        shoebox.data[0, 0, i] = 10.0
        shoebox.background[0, 0, i] = 1.0
    shoebox.mask[0, 0, 0] = MaskCode.Valid | MaskCode.Foreground
    shoebox.mask[0, 0, 1] = MaskCode.Valid | MaskCode.Foreground
    shoebox.mask[0, 0, 2] = (
        MaskCode.Valid | MaskCode.Background | MaskCode.BackgroundUsed
    )
    shoebox.mask[0, 0, 3] = (
        MaskCode.Valid | MaskCode.Background | MaskCode.BackgroundUsed
    )
    shoebox.mask[0, 0, 4] = MaskCode.Valid | MaskCode.Background
    shoebox.mask[0, 0, 5] = MaskCode.Valid | MaskCode.Background

    sliced = FrameSlicedShoebox(
        shoebox, (1, 0, 0), **frame_sliced_shoebox_test_models(1, 1)
    )

    # Two foreground pixels of 10 - 1 each
    assert list(sliced.summation_intensity) == [18.0]
    # Only the two BackgroundUsed pixels count, so m/n is 2/2 and the variance
    # is the sum of the raw foreground values plus the summed background times
    # that ratio, 20 + 2 * 1
    assert list(sliced.summation_intensity_variance) == [22.0]
    assert list(sliced.summation_intensity_valid) == [True]


def test_frame_sliced_shoebox_overlapped_foreground_is_invalid():
    from dials.algorithms.shoebox import MaskCode

    shoebox = Shoebox(0, (0, 2, 0, 1, 0, 2))
    shoebox.allocate_data()
    shoebox.allocate_background()
    for k in range(2):
        for i in range(2):
            shoebox.data[k, 0, i] = 10.0
            shoebox.background[k, 0, i] = 1.0
            shoebox.mask[k, 0, i] = MaskCode.Valid | MaskCode.Foreground

    # A valid foreground pixel on frame 1 that overlaps another reflection
    shoebox.mask[1, 0, 1] |= MaskCode.Overlapped

    sliced = FrameSlicedShoebox(
        shoebox, (1, 0, 0), **frame_sliced_shoebox_test_models(1, 2)
    )

    # An overlapped pixel is excluded from the sum and invalidates its frame,
    # even though it is both valid and foreground
    assert list(sliced.summation_intensity) == [18.0, 9.0]
    assert list(sliced.summation_intensity_valid) == [True, False]


def frame_sliced_shoebox_test_data_with_background_used():
    """The test shoebox, with a row of its background pixels marked as used to
    calculate the background, so that the m/n term of the variance is non-zero"""
    from dials.algorithms.shoebox import MaskCode

    shoebox = frame_sliced_shoebox_test_data()
    for k in range(3):
        for i in range(4):
            shoebox.mask[k, 2, i] |= MaskCode.BackgroundUsed
    return shoebox


def test_frame_sliced_shoebox_summation_matches_summation_integration():
    """The intensity of each frame should equal integrating that frame on its
    own. Its variance should not: see the test below"""
    shoebox = frame_sliced_shoebox_test_data_with_background_used()
    sliced = FrameSlicedShoebox(
        shoebox, (1, 0, 0), **frame_sliced_shoebox_test_models()
    )

    for k in range(3):
        # Build a shoebox holding frame k alone and integrate it
        single = Shoebox(0, (0, 4, 0, 3, 10 + k, 11 + k))
        single.allocate_data()
        single.allocate_background()
        for j in range(3):
            for i in range(4):
                single.data[0, j, i] = shoebox.data[k, j, i]
                single.background[0, j, i] = shoebox.background[k, j, i]
                single.mask[0, j, i] = shoebox.mask[k, j, i]
        expected = single.summed_intensity()

        assert sliced.summation_intensity[k] == pytest.approx(expected.observed.value)
        assert sliced.summation_intensity_valid[k] == expected.observed.success


def test_frame_sliced_shoebox_summation_variances_sum():
    """The background of a frame was determined from the background pixels of
    the whole shoebox, so the m/n term of the variance of a frame uses the pixel
    counts of the whole shoebox rather than those of the frame alone. That makes
    the variance a sum over pixels, so the variances of the frames sum to the
    variance of the whole reflection"""
    shoebox = frame_sliced_shoebox_test_data_with_background_used()
    sliced = FrameSlicedShoebox(
        shoebox, (1, 0, 0), **frame_sliced_shoebox_test_models()
    )

    # Five valid foreground pixels and eleven pixels used for the background, so
    # m/n is 5/11. The foreground values sum to 8 and their background to 2.5
    expected = [2.0 + 1.0 * 5 / 11, 6.0 + 1.5 * 5 / 11, 0.0]
    assert list(sliced.summation_intensity_variance) == pytest.approx(expected)

    # That is the variance of the summation integration of the whole shoebox
    assert sum(sliced.summation_intensity_variance) == pytest.approx(
        shoebox.summed_intensity().observed.variance
    )


def test_frame_sliced_shoebox_excitation_error():
    """The excitation error is the distance of the reciprocal lattice point from
    the surface of the Ewald sphere along the beam, positive if the point is
    inside the sphere"""
    shoebox = frame_sliced_shoebox_test_data()
    models = frame_sliced_shoebox_test_models()

    # The 1 Angstrom beam along -z gives an Ewald sphere of unit radius centred
    # on (0, 0, 1), and the reciprocal lattice point of h is at UB * h

    # (6, 0, 2) is at (0.6, 0, 0.2), which lies exactly on the sphere
    sliced = FrameSlicedShoebox(shoebox, (6, 0, 2), **models)
    assert list(sliced.excitation_error) == pytest.approx([0.0, 0.0, 0.0])

    # (0, 0, 10) is at the centre of the sphere. Moving it a whole radius along
    # the beam, which is -z, brings it to the origin and onto the surface
    sliced = FrameSlicedShoebox(shoebox, (0, 0, 10), **models)
    assert list(sliced.excitation_error) == pytest.approx([1.0, 1.0, 1.0])

    # (1, 0, 0) is at (0.1, 0, 0), which is just outside the sphere
    sliced = FrameSlicedShoebox(shoebox, (1, 0, 0), **models)
    expected = math.sqrt(0.99) - 1.0
    assert expected < 0.0
    assert list(sliced.excitation_error) == pytest.approx([expected] * 3)


def test_frame_sliced_shoebox_excitation_error_is_measured_along_the_beam():
    """The distance is measured along the beam, not radially from the centre of
    the Ewald sphere, so the two differ away from the beam axis"""
    shoebox = frame_sliced_shoebox_test_data()

    # (6, 0, 3) is at (0.6, 0, 0.3), which is inside the unit sphere centred on
    # (0, 0, 1). Moving it 0.1 along the beam brings it to (0.6, 0, 0.2), on the
    # surface, whereas it is only sqrt(0.85) from the centre of the sphere
    sliced = FrameSlicedShoebox(
        shoebox, (6, 0, 3), **frame_sliced_shoebox_test_models()
    )
    assert list(sliced.excitation_error) == pytest.approx([0.1] * 3)
    assert 1.0 - math.sqrt(0.85) != pytest.approx(0.1)


def test_frame_sliced_shoebox_excitation_error_can_be_undefined():
    """A point further off the beam axis than the radius of the Ewald sphere is
    not brought to the surface by any movement along the beam"""
    shoebox = frame_sliced_shoebox_test_data()

    # (12, 0, 5) is 1.2 from the beam axis, which the unit sphere never reaches
    sliced = FrameSlicedShoebox(
        shoebox, (12, 0, 5), **frame_sliced_shoebox_test_models()
    )
    assert all(math.isnan(e) for e in sliced.excitation_error)


def test_frame_sliced_shoebox_excitation_error_is_calculated_per_frame():
    """Each frame uses the orientation matrix and beam vector of its own image"""
    shoebox = frame_sliced_shoebox_test_data()
    models = frame_sliced_shoebox_test_models()

    # A crystal cell and a wavelength that both differ on each of the frames
    cells = (10.0, 11.0, 12.0)
    wavelengths = (1.0, 1.1, 1.2)
    models["UB"] = flex.mat3_double(
        [(1 / a, 0, 0, 0, 1 / a, 0, 0, 0, 1 / a) for a in cells]
    )
    models["s0"] = flex.vec3_double([(0, 0, -1 / w) for w in wavelengths])

    sliced = FrameSlicedShoebox(shoebox, (1, 0, 0), **models)

    expected = [
        math.sqrt((1 / w) ** 2 - (1 / a) ** 2) - 1 / w
        for a, w in zip(cells, wavelengths)
    ]
    assert list(sliced.excitation_error) == pytest.approx(expected)


def test_frame_sliced_shoebox_partiality():
    """The partiality of a frame is the portion of the rocking curve of the
    reflection that lies between the boundaries of that frame"""
    shoebox = frame_sliced_shoebox_test_data()

    # A rocking curve centred on the boundary between the first two frames of
    # the shoebox, which are at 1.0 -> 1.1 -> 1.2 -> 1.3 radians
    sliced = FrameSlicedShoebox(
        shoebox, (1, 0, 0), **frame_sliced_shoebox_test_models(phi_cal=1.1)
    )
    assert list(sliced.partiality) == pytest.approx(
        [
            gaussian_cdf(0.0) - gaussian_cdf(-1.0),
            gaussian_cdf(1.0) - gaussian_cdf(0.0),
            gaussian_cdf(2.0) - gaussian_cdf(1.0),
        ]
    )

    # The rocking curve is symmetric about the boundary it is centred on, so the
    # two frames either side of it take equal shares of the reflection. That
    # holds for the boundaries between frames, not for the frame centres
    assert sliced.partiality[0] == pytest.approx(sliced.partiality[1])
    assert sliced.partiality[1] != pytest.approx(sliced.partiality[2])

    # A rocking curve much narrower than a frame is recorded on one frame alone
    sliced = FrameSlicedShoebox(
        shoebox,
        (1, 0, 0),
        **frame_sliced_shoebox_test_models(phi_cal=1.15, sigma_phi=0.001),
    )
    assert list(sliced.partiality) == pytest.approx([0.0, 1.0, 0.0], abs=1e-12)


def test_frame_sliced_shoebox_partiality_sums_over_the_frames():
    """The partialities of the frames sum to the partiality of the whole
    shoebox, which is what PartialityCalculator3D returns"""
    shoebox = frame_sliced_shoebox_test_data()

    for phi_cal, sigma_phi in ((1.15, 0.1), (1.0, 0.05), (0.8, 0.3), (2.0, 0.2)):
        sliced = FrameSlicedShoebox(
            shoebox,
            (1, 0, 0),
            **frame_sliced_shoebox_test_models(phi_cal=phi_cal, sigma_phi=sigma_phi),
        )

        # The shoebox spans images 11 to 13, that is 1.0 to 1.3 radians
        whole = gaussian_cdf((1.3 - phi_cal) / sigma_phi) - gaussian_cdf(
            (1.0 - phi_cal) / sigma_phi
        )
        assert sum(sliced.partiality) == pytest.approx(whole)
        assert all(0.0 <= p <= 1.0 for p in sliced.partiality)


def test_frame_sliced_shoebox_partiality_needs_a_rocking_curve_width():
    shoebox = frame_sliced_shoebox_test_data()

    for sigma_phi in (0.0, -0.1):
        with pytest.raises(RuntimeError):
            FrameSlicedShoebox(
                shoebox,
                (1, 0, 0),
                **frame_sliced_shoebox_test_models(sigma_phi=sigma_phi),
            )


def test_frame_sliced_shoebox_models_are_looked_up_by_image_number():
    """The model arrays cover a whole scan, not just the frames of the shoebox"""
    shoebox = frame_sliced_shoebox_test_data()

    # Models for the whole of a 20 image scan, of which the shoebox covers 11-13,
    # with a crystal cell that shrinks as the scan goes on
    models = frame_sliced_shoebox_test_models(1, 20)
    cells = [100.0 / f for f in range(1, 21)]
    models["UB"] = flex.mat3_double(
        [(1 / a, 0, 0, 0, 1 / a, 0, 0, 0, 1 / a) for a in cells]
    )

    sliced = FrameSlicedShoebox(shoebox, (1, 0, 0), **models)
    assert list(sliced.phi) == pytest.approx([1.05, 1.15, 1.25])
    expected = [math.sqrt(1.0 - (f / 100.0) ** 2) - 1.0 for f in (11, 12, 13)]
    assert list(sliced.excitation_error) == pytest.approx(expected)

    # The scan does not reach as far as the last frame of the shoebox
    with pytest.raises(RuntimeError):
        FrameSlicedShoebox(
            shoebox, (1, 0, 0), **frame_sliced_shoebox_test_models(1, 12)
        )

    # The scan starts after the first frame of the shoebox
    with pytest.raises(RuntimeError):
        FrameSlicedShoebox(
            shoebox, (1, 0, 0), **frame_sliced_shoebox_test_models(12, 3)
        )


def test_frame_sliced_shoebox_requires_consistent_model_arrays():
    """Every model array must cover the same images"""
    shoebox = frame_sliced_shoebox_test_data()

    for name in ("UB", "s0", "phi_scan_points"):
        models = frame_sliced_shoebox_test_models(1, 20)
        models[name] = models[name][:-1]
        with pytest.raises(RuntimeError):
            FrameSlicedShoebox(shoebox, (1, 0, 0), **models)

    # There must be one more scan point than there are frames
    models = frame_sliced_shoebox_test_models(1, 20)
    models["phi"] = models["phi_scan_points"]
    with pytest.raises(RuntimeError):
        FrameSlicedShoebox(shoebox, (1, 0, 0), **models)


def test_frame_sliced_shoebox_arrays_are_read_only():
    sliced = FrameSlicedShoebox(
        frame_sliced_shoebox_test_data(),
        (1, 0, 0),
        **frame_sliced_shoebox_test_models(),
    )

    for name in (
        "frames",
        "phi",
        "excitation_error",
        "partiality",
        "summation_intensity",
        "summation_intensity_variance",
        "summation_intensity_valid",
    ):
        with pytest.raises(AttributeError):
            setattr(sliced, name, None)


def test_frame_sliced_shoebox_requires_allocated_arrays():
    models = frame_sliced_shoebox_test_models()

    # Neither the data nor the background are allocated
    shoebox = Shoebox(0, (0, 4, 0, 3, 10, 13))
    with pytest.raises(RuntimeError):
        FrameSlicedShoebox(shoebox, (1, 0, 0), **models)

    # The data is allocated, but the background is not
    shoebox.allocate_data()
    with pytest.raises(RuntimeError):
        FrameSlicedShoebox(shoebox, (1, 0, 0), **models)

    shoebox.allocate_background()
    assert len(FrameSlicedShoebox(shoebox, (1, 0, 0), **models)) == 3


def test_frame_sliced_shoebox_rejects_flat_shoebox():
    shoebox = frame_sliced_shoebox_test_data()
    shoebox.flatten()
    with pytest.raises(RuntimeError):
        FrameSlicedShoebox(shoebox, (1, 0, 0), **frame_sliced_shoebox_test_models())

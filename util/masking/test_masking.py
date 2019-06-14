from __future__ import absolute_import, division, print_function

import math
import pytest

from scitbx.array_family import flex

from dxtbx.model.goniometer import GoniometerFactory
from dxtbx.model.detector import DetectorFactory
from dials.util.masking import GoniometerShadowMaskGenerator


@pytest.fixture
def kappa_goniometer():
    def _construct_kappa_goniometer(phi, kappa, omega):
        phi_axis = (1.0, 0.0, 0.0)
        kappa_axis = (0.914, 0.279, -0.297)
        omega_axis = (1.0, 0.0, 0.0)
        axes = flex.vec3_double((phi_axis, kappa_axis, omega_axis))
        angles = flex.double((phi, kappa, omega))
        names = flex.std_string(("GON_PHI", "GON_KAPPA", "GON_OMEGA"))
        return GoniometerFactory.make_multi_axis_goniometer(
            axes, angles, names, scan_axis=2
        )

    return _construct_kappa_goniometer


@pytest.fixture
def pilatus_6M():
    def _construct_detector(distance):
        return DetectorFactory.simple(
            sensor="PAD",
            distance=distance,
            beam_centre=(216.87, 211.32),
            fast_direction="+x",
            slow_direction="-y",
            pixel_size=(0.172, 0.172),
            image_size=(2463, 2527),
            trusted_range=(-1, 1e8),
        )

    return _construct_detector


def test_GoniometerShadowMaskGenerator(kappa_goniometer, pilatus_6M):
    goniometer = kappa_goniometer(phi=0, kappa=180, omega=0)
    # Simple model of cone around goniometer phi axis
    # Exact values don't matter, only the ratio of height/radius
    height = 50  # mm
    radius = 20  # mm

    steps_per_degree = 1
    theta = (
        flex.double([range(360 * steps_per_degree)])
        * math.pi
        / 180
        * 1
        / steps_per_degree
    )
    y = radius * flex.cos(theta)
    z = radius * flex.sin(theta)
    x = flex.double(theta.size(), height)

    coords = flex.vec3_double(zip(x, y, z))
    coords.insert(0, (0, 0, 0))

    masker = GoniometerShadowMaskGenerator(
        goniometer, coords, flex.size_t(len(coords), 0)
    )

    detector = pilatus_6M(distance=170)
    shadow = masker.project_extrema(detector, scan_angle=0)
    assert len(shadow) == len(detector)
    assert len(shadow[0]) == 145
    assert shadow[0][0] == pytest.approx((1752.4703959222845, 672.2428154208985))
    extrema = masker.extrema_at_scan_angle(scan_angle=0)
    assert len(extrema) == len(coords)
    assert extrema[1] == pytest.approx(
        (43.60448791048147, 8.572923552543028, -30.416337975287743)
    )
    mask = masker.get_mask(detector, 0)
    assert len(mask) == len(detector)
    assert mask[0].all() == tuple(reversed(detector[0].get_image_size()))
    assert mask[0].count(True) == 5570865

    shadow = masker.project_extrema(detector, scan_angle=-45)
    assert len(shadow[0]) == 226
    assert shadow[0][0] == pytest.approx((1623.5133257425446, 1254.0966982780835))
    extrema = masker.extrema_at_scan_angle(scan_angle=-45)
    assert len(extrema) == len(coords)
    assert extrema[1] == pytest.approx(
        (43.60448791048147, -15.445626462590827, -27.569571219784905)
    )
    mask = masker.get_mask(detector, -45)
    assert mask[0].all() == tuple(reversed(detector[0].get_image_size()))
    assert mask[0].count(True) == 5467810

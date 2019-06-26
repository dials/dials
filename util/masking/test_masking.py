from __future__ import absolute_import, division, print_function

import math
import pytest

from scitbx.array_family import flex

from dxtbx.model.goniometer import GoniometerFactory
from dxtbx.model.detector import DetectorFactory
from dials.util.masking import (
    GoniometerShadowMaskGenerator,
    PyGoniometerShadowMaskGenerator,
)


@pytest.fixture
def kappa_goniometer():
    def _construct_goniometer(phi, kappa, omega):
        phi_axis = (1.0, 0.0, 0.0)
        kappa_axis = (0.914, 0.279, -0.297)
        omega_axis = (1.0, 0.0, 0.0)
        axes = flex.vec3_double((phi_axis, kappa_axis, omega_axis))
        angles = flex.double((phi, kappa, omega))
        names = flex.std_string(("GON_PHI", "GON_KAPPA", "GON_OMEGA"))
        return GoniometerFactory.make_multi_axis_goniometer(
            axes, angles, names, scan_axis=2
        )

    return _construct_goniometer


@pytest.fixture
def smargon_goniometer():
    def _construct_goniometer(phi, chi, omega):
        phi_axis = (1.0, 0.0, 0.0)
        chi_axis = (0, 0, -1)
        omega_axis = (1.0, 0.0, 0.0)
        axes = flex.vec3_double((phi_axis, chi_axis, omega_axis))
        angles = flex.double((phi, chi, omega))
        names = flex.std_string(("GON_PHI", "GON_CHI", "GON_OMEGA"))
        return GoniometerFactory.make_multi_axis_goniometer(
            axes, angles, names, scan_axis=2
        )

    return _construct_goniometer


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


@pytest.fixture(params=["cpp", "python"])
def kappa_goniometer_shadow_masker(request):
    def _construct_shadow_masker(goniometer):
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
        y = radius * flex.cos(-theta)
        z = radius * flex.sin(-theta)
        x = flex.double(theta.size(), height)

        coords = flex.vec3_double(zip(x, y, z))
        coords.append(coords[0])  # ensure closure of polygon
        coords.insert(0, (0, 0, 0))

        # from matplotlib import pyplot
        # pyplot.plot(y, z)
        # pyplot.show()

        if request.param == "python":
            return PyGoniometerShadowMaskGenerator(
                goniometer, coords, flex.size_t(len(coords), 0)
            )
        return GoniometerShadowMaskGenerator(
            goniometer, coords, flex.size_t(len(coords), 0)
        )

    return _construct_shadow_masker


@pytest.fixture(params=["cpp", "python"])
def smargon_shadow_masker(request):
    def _construct_shadow_masker(goniometer):
        from dxtbx.format.SmarGonShadowMask import (
            SmarGonShadowMaskGenerator,
            PySmarGonShadowMaskGenerator,
        )

        if request.param == "python":
            return PySmarGonShadowMaskGenerator(goniometer)
        return SmarGonShadowMaskGenerator(goniometer)

    return _construct_shadow_masker


def test_GoniometerShadowMaskGenerator_kappa_180_omega_0(
    kappa_goniometer, pilatus_6M, kappa_goniometer_shadow_masker
):
    goniometer = kappa_goniometer(phi=0, kappa=180, omega=0)
    detector = pilatus_6M(distance=170)
    masker = kappa_goniometer_shadow_masker(goniometer)

    scan_angle = 0  # omega = 0, shadow in top right corner of detector
    extrema = masker.extrema_at_scan_angle(scan_angle)
    assert len(extrema) == 362
    assert extrema[1] == pytest.approx(
        (43.60448791048147, 8.572923552543028, -30.416337975287743)
    )

    shadow = masker.project_extrema(detector, scan_angle)
    assert len(shadow) == len(detector)
    assert len(shadow[0]) == 145

    mask = masker.get_mask(detector, scan_angle)
    assert len(mask) == len(detector)
    assert mask[0].all() == tuple(reversed(detector[0].get_image_size()))
    assert mask[0].count(True) == pytest.approx(5570865)


def test_GoniometerShadowMaskGenerator_kappa_180_omega_m45(
    kappa_goniometer, pilatus_6M, kappa_goniometer_shadow_masker
):
    goniometer = kappa_goniometer(phi=0, kappa=180, omega=0)
    detector = pilatus_6M(distance=170)
    masker = kappa_goniometer_shadow_masker(goniometer)

    scan_angle = -45  # omega = -45, shadow on centre-right of detector
    extrema = masker.extrema_at_scan_angle(scan_angle)
    assert len(extrema) == 362
    assert extrema[1] == pytest.approx(
        (43.60448791048147, -15.445626462590827, -27.569571219784905)
    )
    shadow = masker.project_extrema(detector, scan_angle)
    assert len(shadow[0]) in (226, 227)
    # assert shadow[0][0] == pytest.approx((1623.5133257425446, 1254.0966982780835))
    mask = masker.get_mask(detector, scan_angle)
    assert mask[0].all() == tuple(reversed(detector[0].get_image_size()))
    assert mask[0].count(True) == pytest.approx(5467810)


def test_GoniometerShadowMaskGenerator_kappa_180_omega_p45(
    kappa_goniometer, pilatus_6M, kappa_goniometer_shadow_masker
):
    goniometer = kappa_goniometer(phi=0, kappa=180, omega=0)
    detector = pilatus_6M(distance=170)
    masker = kappa_goniometer_shadow_masker(goniometer)

    scan_angle = 45  # omega = +45, no shadow on detector
    shadow = masker.project_extrema(detector, scan_angle)
    assert len(shadow[0]) == 0
    extrema = masker.extrema_at_scan_angle(scan_angle)
    assert len(extrema) == 362
    assert extrema[1] == pytest.approx(
        (43.60448791048147, 27.56957121978491, -15.445626462590822)
    )
    mask = masker.get_mask(detector, scan_angle)
    assert mask[0] is None or mask[0].count(False) == 0


def test_GoniometerShadowMaskGenerator_kappa_m70_omega_p100(
    kappa_goniometer, pilatus_6M, kappa_goniometer_shadow_masker
):
    # goniometer shadow does not intersect with detector panel
    goniometer = kappa_goniometer(phi=0, kappa=-70, omega=0)
    detector = pilatus_6M(distance=170)
    masker = kappa_goniometer_shadow_masker(goniometer)

    scan_angle = 100
    extrema = masker.extrema_at_scan_angle(scan_angle)
    assert len(extrema) == 362
    assert extrema[1] == pytest.approx(
        (42.318198019878935, 8.61724085750561, 32.17006801910824)
    )
    shadow = masker.project_extrema(detector, scan_angle)
    assert len(shadow[0]) == 0
    mask = masker.get_mask(detector, scan_angle)
    assert mask[0] is None or mask[0].count(False) == 0


def test_SmarGonShadowMaskGenerator_p48_c45_o95(
    smargon_goniometer, pilatus_6M, smargon_shadow_masker
):
    goniometer = smargon_goniometer(phi=48, chi=45, omega=100)
    detector = pilatus_6M(distance=170)
    masker = smargon_shadow_masker(goniometer)

    scan_angle = 100
    extrema = masker.extrema_at_scan_angle(scan_angle)
    assert len(extrema) == 82
    assert extrema[1] == pytest.approx(
        (22.106645739466337, 13.963679418895005, -22.47914302717216)
    )
    shadow = masker.project_extrema(detector, scan_angle)
    assert len(shadow[0]) in (15, 16)
    mask = masker.get_mask(detector, scan_angle)
    assert mask[0].count(True) == pytest.approx(5716721, 3e-5)


def test_SmarGonShadowMaskGenerator_p0_c90_o50(
    smargon_goniometer, pilatus_6M, smargon_shadow_masker
):
    for phi in (0, 90, 180):
        # There should be no dependence of the shadow on phi
        goniometer = smargon_goniometer(phi=phi, chi=90, omega=50)
        detector = pilatus_6M(distance=214.75)
        masker = smargon_shadow_masker(goniometer)

        scan_angle = 50
        extrema = masker.extrema_at_scan_angle(scan_angle)
        assert len(extrema) == 82
        assert extrema[1] == pytest.approx(
            (-1.7364817766692957, -13.66792605230091, -31.609688838521162)
        )
        shadow = masker.project_extrema(detector, scan_angle)
        assert len(shadow[0]) == 11
        mask = masker.get_mask(detector, scan_angle)
        assert mask[0].count(True) == pytest.approx(4614865, 2e-4)

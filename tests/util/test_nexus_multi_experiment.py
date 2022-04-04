from __future__ import annotations

import math
import os
import random

import pytest

from dxtbx.model.experiment_list import ExperimentListFactory
from scitbx import matrix

import dials.util.nexus
from dials.command_line.export import phil_scope
from dials.util.nexus.nx_mx import (
    polarization_normal_to_stokes,
    polarization_stokes_to_normal,
)


def test_polarization_conversion():
    def conv(n, p):
        S0, S1, S2, S3 = polarization_normal_to_stokes(n, p)
        return polarization_stokes_to_normal(S0, S1, S2, S3)

    for i in range(1000):
        n1 = matrix.col((random.uniform(-1, 1), random.uniform(-1, 1), 0)).normalize()
        p1 = random.uniform(0, 1)
        n2, p2 = conv(n1, p1)
        angle = n2.angle(n1)
        assert (
            abs(angle) == pytest.approx(0, abs=1e-7)
            or angle == pytest.approx(math.pi)
            or angle == pytest.approx(math.pi * 2)
        )
        assert p1 == pytest.approx(p2, abs=1e-7)


def run_single(experiments1, filename):
    # Dump the file
    params = phil_scope.extract().nxs
    params.hklout = filename
    dials.util.nexus.dump(experiments1, None, params)

    # Load the file
    experiments2, reflections = dials.util.nexus.load(filename)
    assert experiments2 is not None
    assert reflections is None
    assert len(experiments2) == len(experiments1)

    index1 = []
    index2 = []

    for exp1, exp2 in zip(experiments1, experiments2):

        # Check the beam
        b1 = exp1.beam
        b2 = exp2.beam
        assert all(
            d1 == pytest.approx(d2)
            for d1, d2 in zip(
                b1.get_sample_to_source_direction(), b2.get_sample_to_source_direction()
            )
        )
        assert b1.get_wavelength() == pytest.approx(b2.get_wavelength())
        assert b1.get_polarization_fraction() == pytest.approx(
            b2.get_polarization_fraction()
        )
        n1 = matrix.col(b1.get_polarization_normal())
        n2 = matrix.col(b2.get_polarization_normal())
        angle = n2.angle(n1)
        assert (
            angle == pytest.approx(0, abs=1e-7)
            or angle == pytest.approx(math.pi)
            or angle == pytest.approx(2 * math.pi)
        )

        # Check the goniometer
        g1 = exp1.goniometer
        g2 = exp2.goniometer
        if g1 is None:
            assert g2 is None
        else:
            assert all(
                d1 == pytest.approx(d2)
                for d1, d2 in zip(g1.get_rotation_axis(), g2.get_rotation_axis())
            )
            assert all(
                d1 == pytest.approx(d2)
                for d1, d2 in zip(g1.get_fixed_rotation(), g2.get_fixed_rotation())
            )
            assert all(
                d1 == pytest.approx(d2)
                for d1, d2 in zip(g1.get_setting_rotation(), g2.get_setting_rotation())
            )

        # Check the scan
        s1 = exp1.scan
        s2 = exp2.scan
        if s1 is None:
            assert s2 is None
        else:
            assert len(s1) == len(s2)
            assert s1.get_image_range() == s2.get_image_range()
            assert s1.get_oscillation() == pytest.approx(s2.get_oscillation())
            assert s1.get_exposure_times() == pytest.approx(s2.get_exposure_times())
            assert s1.get_epochs() == pytest.approx(s2.get_epochs())

        # Check the detector
        d1 = exp1.detector
        d2 = exp2.detector
        assert len(d1) == len(d2)
        for p1, p2 in zip(d1, d2):
            assert p1.get_type() == p2.get_type()
            assert p1.get_material() == p2.get_material()
            assert p1.get_thickness() == p2.get_thickness()
            assert p1.get_image_size() == p2.get_image_size()
            assert p1.get_pixel_size() == p2.get_pixel_size()
            assert p1.get_trusted_range() == p2.get_trusted_range()
            assert p1.get_fast_axis() == pytest.approx(p2.get_fast_axis())
            assert p1.get_slow_axis() == pytest.approx(p2.get_slow_axis())
            assert p1.get_origin() == pytest.approx(p2.get_origin())

        # Check the crystal
        c1 = exp1.crystal
        c2 = exp2.crystal
        assert c1.get_space_group() == c2.get_space_group()
        assert c1.get_unit_cell().parameters() == pytest.approx(
            c2.get_unit_cell().parameters()
        )
        assert c1.get_A() == pytest.approx(c2.get_A())
        assert c1.num_scan_points == c2.num_scan_points
        for i in range(c1.num_scan_points):
            assert c1.get_A_at_scan_point(i) == pytest.approx(c2.get_A_at_scan_point(i))
            uc1 = c1.get_unit_cell_at_scan_point(i)
            uc2 = c2.get_unit_cell_at_scan_point(i)
            assert uc1.parameters() == pytest.approx(uc2.parameters())

        index1.append(
            (
                id(exp1.beam),
                id(exp1.detector),
                id(exp1.goniometer),
                id(exp1.scan),
                id(exp1.crystal),
            )
        )

        index2.append(
            (
                id(exp2.beam),
                id(exp2.detector),
                id(exp2.goniometer),
                id(exp2.scan),
                id(exp2.crystal),
            )
        )

    # Get a list of all beam etc
    beam1, detector1, goniometer1, scan1, crystal1 = zip(*index1)
    beam2, detector2, goniometer2, scan2, crystal2 = zip(*index2)

    # If any models are shared then check they are shared in both
    num = len(beam1)
    for i in range(0, num - 1):
        for j in range(1, num):
            assert (beam1[i] == beam1[j]) is (beam2[i] == beam2[j])
            assert (detector1[i] == detector1[j]) is (detector2[i] == detector2[j])
            assert (goniometer1[i] == goniometer1[j]) is (
                goniometer2[i] == goniometer2[j]
            )
            assert (scan1[i] == scan1[j]) is (scan2[i] == scan2[j])
            assert (crystal1[i] == crystal1[j]) is (crystal2[i] == crystal2[j])


@pytest.mark.parametrize(
    "filename",
    [
        "single",
        "multiple_unrelated",
        "multi_crystal",
        "two_colour",
        "multiple_sweeps",
        "stills",
    ],
)
def test_nexus_dump_and_reload(dials_regression, tmpdir, filename):
    path = os.path.join(dials_regression, "nexus_test_data", "shared_models")
    filename_in = os.path.join(path, f"{filename}.json")
    filename_out = tmpdir.join(f"{filename}.nxs").strpath
    experiments = ExperimentListFactory.from_json_file(filename_in)
    run_single(experiments, filename_out)

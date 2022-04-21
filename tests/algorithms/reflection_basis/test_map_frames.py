from __future__ import annotations

import math
import random

from dxtbx.serialize import load
from scitbx import matrix
from scitbx.array_family import flex

from dials.algorithms.profile_model.gaussian_rs import (
    BBoxCalculator3D,
    CoordinateSystem,
)
from dials.algorithms.profile_model.gaussian_rs.transform import (
    MapFramesForward,
    MapFramesReverse,
)


def test_map_frames_forward(dials_data):
    sequence = load.imageset(
        dials_data("centroid_test_data", pathlib=True) / "sweep.json"
    )

    # Get the models
    beam = sequence.get_beam()
    detector = sequence.get_detector()
    gonio = sequence.get_goniometer()
    scan = sequence.get_scan()

    # Set the delta_divergence/mosaicity
    n_sigma = 3
    sigma_divergence = 0.060 * math.pi / 180
    mosaicity = 0.154 * math.pi / 180
    delta_divergence = n_sigma * sigma_divergence
    delta_mosaicity = n_sigma * mosaicity

    # Set the grid size
    grid_size = (4, 4, 4)

    # Create the E3 fraction object
    transform = MapFramesForward(
        scan.get_array_range()[0],
        scan.get_oscillation(deg=False)[0],
        scan.get_oscillation(deg=False)[1],
        mosaicity,
        n_sigma,
        grid_size[2],
    )

    # Create the bounding box calculator
    calculate_bbox = BBoxCalculator3D(
        beam, detector, gonio, scan, delta_divergence, delta_mosaicity
    )

    assert len(detector) == 1
    s0 = beam.get_s0()
    m2 = gonio.get_rotation_axis()
    s0_length = matrix.col(beam.get_s0()).length()

    for i in range(100):

        # Get random x, y, z
        x = random.uniform(0, 2000)
        y = random.uniform(0, 2000)
        z = random.uniform(0, 9)

        # Get random s1, phi, panel
        s1 = matrix.col(detector[0].get_pixel_lab_coord((x, y))).normalize() * s0_length
        phi = scan.get_angle_from_array_index(z, deg=False)
        panel = 0

        # Calculate the bounding box
        bbox = calculate_bbox(s1, z, panel)

        # Create the XDS coordinate system
        xcs = CoordinateSystem(m2, s0, s1, phi)

        # Calculate the transform fraction
        fraction = transform(bbox[4:], phi, xcs.zeta())

        # Ensure the minimum and maximum are 0 < 1
        fmax = flex.max(fraction)
        fmin = flex.min(fraction)
        assert fmax <= (1.0 + 5e-15) and fmax > 0.0, f"{fmax:.16f} not between 0 and 1"
        assert fmin >= 0.0 and fmin <= 1.0

        # Ensure the fraction for each image frame adds up to 1.0 for
        # all those frames completely within the grid
        for j in range(1, fraction.all()[0] - 1):
            tot = flex.sum(fraction[j : j + 1, :])
            assert abs(tot - 1.0) < 1e-7

        # Ensure the frames follow a progression through the grid. I.e,
        # check that values increase then decrease and don't jump around
        for j in range(fraction.all()[0]):
            f = fraction[j : j + 1, :]
            last = f[0]
            rev = False
            for i in range(1, len(f)):
                curr = f[1]
                if rev is False:
                    if curr < last:
                        rev = True
                else:
                    assert curr <= last
                last = curr


def test_map_frames_reverse(dials_data):
    sequence = load.imageset(
        dials_data("centroid_test_data", pathlib=True) / "sweep.json"
    )

    # Get the models
    beam = sequence.get_beam()
    detector = sequence.get_detector()
    gonio = sequence.get_goniometer()
    scan = sequence.get_scan()

    # Set the delta_divergence/mosaicity
    n_sigma = 3
    sigma_divergence = 0.060 * math.pi / 180
    mosaicity = 0.154 * math.pi / 180
    delta_divergence = n_sigma * sigma_divergence
    delta_mosaicity = n_sigma * mosaicity

    # Set the grid size
    grid_size = (4, 4, 4)

    # Create the E3 fraction object
    transform = MapFramesReverse(
        scan.get_array_range()[0],
        scan.get_oscillation(deg=False)[0],
        scan.get_oscillation(deg=False)[1],
        mosaicity,
        n_sigma,
        grid_size[2],
    )

    # Create the bounding box calculator
    calculate_bbox = BBoxCalculator3D(
        beam, detector, gonio, scan, delta_divergence, delta_mosaicity
    )

    s0 = beam.get_s0()
    m2 = gonio.get_rotation_axis()
    s0_length = matrix.col(beam.get_s0()).length()

    for i in range(100):

        # Get random x, y, z
        x = random.uniform(0, 2000)
        y = random.uniform(0, 2000)
        z = random.uniform(0, 9)

        # Get random s1, phi, panel
        s1 = matrix.col(detector[0].get_pixel_lab_coord((x, y))).normalize() * s0_length
        phi = scan.get_angle_from_array_index(z, deg=False)
        panel = 0

        # Calculate the bounding box
        bbox = calculate_bbox(s1, phi, panel)
        x1, x2 = bbox[0], bbox[1]
        y1, y2 = bbox[2], bbox[3]
        z1, z2 = bbox[4], bbox[5]
        if x1 == 0 or y1 == 0 or z1 == 0:
            continue
        if x2 == 2000 or y2 == 2000 or z2 == 9:
            continue

        # Create the XDS coordinate system
        xcs = CoordinateSystem(m2, s0, s1, phi)

        # Calculate the transform fraction
        fraction = transform(bbox[4:], phi, xcs.zeta())

        # Ensure the minimum and maximum are 0 < 1
        fmax = flex.max(fraction)
        fmin = flex.min(fraction)
        assert fmax <= 1.0 and fmax > 0.0
        assert fmin >= 0.0 and fmin <= 1.0

        # Ensure the fraction for image adds up to 1.0 for
        # all those images completely within the image
        for v3 in range(fraction.all()[0]):
            tot = flex.sum(fraction[v3 : v3 + 1, :])
            assert abs(tot - 1.0) < 1e-7

        # Ensure the frames follow a progression through the grid. I.e,
        # check that values increase then decrease and don't jump around
        for v3 in range(fraction.all()[0]):
            f = fraction[v3 : v3 + 1, :]
            last = f[0]
            rev = False
            for i in range(1, len(f)):
                curr = f[1]
                if rev is False:
                    if curr < last:
                        rev = True
                else:
                    assert curr <= last
                last = curr


def test_map_forward_reverse(dials_data):
    sequence = load.imageset(
        dials_data("centroid_test_data", pathlib=True) / "sweep.json"
    )

    # Get the models
    beam = sequence.get_beam()
    detector = sequence.get_detector()
    gonio = sequence.get_goniometer()
    scan = sequence.get_scan()

    # Set the delta_divergence/mosaicity
    n_sigma = 3
    sigma_divergence = 0.060 * math.pi / 180
    mosaicity = 0.154 * math.pi / 180
    delta_divergence = n_sigma * sigma_divergence
    delta_mosaicity = n_sigma * mosaicity

    # Set the grid size
    grid_size = (4, 4, 4)

    # Create the E3 fraction object
    transform_forward = MapFramesForward(
        scan.get_array_range()[0],
        scan.get_oscillation(deg=False)[0],
        scan.get_oscillation(deg=False)[1],
        mosaicity,
        n_sigma,
        grid_size[2],
    )

    # Create the E3 fraction object
    transform_reverse = MapFramesReverse(
        scan.get_array_range()[0],
        scan.get_oscillation(deg=False)[0],
        scan.get_oscillation(deg=False)[1],
        mosaicity,
        n_sigma,
        grid_size[2],
    )

    # Create the bounding box calculator
    calculate_bbox = BBoxCalculator3D(
        beam, detector, gonio, scan, delta_divergence, delta_mosaicity
    )

    s0 = beam.get_s0()
    m2 = gonio.get_rotation_axis()
    s0_length = matrix.col(beam.get_s0()).length()

    for i in range(100):

        # Get random x, y, z
        x = random.uniform(0, 2000)
        y = random.uniform(0, 2000)
        z = random.uniform(0, 9)

        # Get random s1, phi, panel
        s1 = matrix.col(detector[0].get_pixel_lab_coord((x, y))).normalize() * s0_length
        phi = scan.get_angle_from_array_index(z, deg=False)
        panel = 0

        # Calculate the bounding box
        bbox = calculate_bbox(s1, phi, panel)

        # Create the XDS coordinate system
        xcs = CoordinateSystem(m2, s0, s1, phi)

        # Calculate the transform fraction
        forward_fraction = transform_forward(bbox[4:], phi, xcs.zeta())

        # Calculate the transform fraction
        reverse_fraction = transform_reverse(bbox[4:], phi, xcs.zeta())

        # Check the same points are non-zero
        eps = 1e-7
        for j in range(forward_fraction.all()[0]):
            for i in range(forward_fraction.all()[1]):
                if forward_fraction[j, i] > 0.0:
                    assert reverse_fraction[i, j] > 0.0
                else:
                    assert reverse_fraction[i, j] < eps

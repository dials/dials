from __future__ import annotations

import math
import pickle
import random
import shutil
import subprocess

import mrcfile
import numpy as np
import pytest
from skimage.measure import EllipseModel

from dxtbx.format.Format import Reader
from dxtbx.imageset import ImageSet, ImageSetData
from dxtbx.model.beam import Beam
from dxtbx.model.detector import Detector
from dxtbx.model.experiment_list import ExperimentList, ExperimentListFactory
from libtbx import easy_run
from scitbx import matrix

from dials.array_family import flex
from dials.command_line.generate_distortion_maps import (
    circle_to_ellipse_transform,
    ellipse_to_circle_transform,
)
from dials.util.image_viewer.slip_viewer.ellipse_frame import extract_ellipse_parameters


def make_detector():
    """Make a dummy 4 panel detector"""
    pixel_size_x = 0.1
    pixel_size_y = 0.1
    npixels_per_panel_x = 120
    npixels_per_panel_y = 160
    distance = 100
    fast = matrix.col((1, 0, 0))
    slow = matrix.col((0, -1, 0))
    s0u = matrix.col((0, 0, -1))

    shift_x = -1.0 * npixels_per_panel_x * pixel_size_x * fast
    shift_y = -1.0 * npixels_per_panel_y * pixel_size_y * slow
    beam_centre = distance * s0u
    orig = beam_centre + shift_x + shift_y

    d = Detector()
    root = d.hierarchy()
    root.set_local_frame(fast.elems, slow.elems, orig.elems)

    local_fast = matrix.col((1.0, 0.0, 0.0))
    local_slow = matrix.col((0.0, 1.0, 0.0))
    span_x = npixels_per_panel_x * pixel_size_x * local_fast
    span_y = npixels_per_panel_y * pixel_size_y * local_slow

    local_origins = [matrix.col((0, 0, 0)), span_x, span_y, span_x + span_y]

    for i, lo in enumerate(local_origins):
        p = d.add_panel()
        p.set_image_size((npixels_per_panel_x, npixels_per_panel_y))
        p.set_pixel_size((pixel_size_x, pixel_size_y))
        p.set_local_frame(local_fast.elems, local_slow.elems, lo.elems)

    return d


def test_translate(dials_data, run_in_tmp_path):
    """Test as written in https://github.com/dials/dials/issues/471. This
    is pretty slow!"""

    data_dir = dials_data("centroid_test_data")
    image_path = data_dir / "centroid_0001.cbf"

    # Generate distortion maps
    cmd = f"dials.generate_distortion_maps {image_path} dx=1 dy=2"
    easy_run.fully_buffered(command=cmd).raise_if_errors()

    # Import without correction
    cmd = f"dials.import {image_path}"
    easy_run.fully_buffered(command=cmd).raise_if_errors()
    expt1 = ExperimentListFactory.from_serialized_format("imported.expt")[0]

    # Should be no dx/dy lookup files
    assert not expt1.imageset.external_lookup.dx.filename
    assert not expt1.imageset.external_lookup.dy.filename

    # Import with correction
    cmd = f"dials.import {image_path} dx=dx.pickle dy=dy.pickle output.experiments=corrected.expt"

    easy_run.fully_buffered(command=cmd).raise_if_errors()
    expt2 = ExperimentListFactory.from_serialized_format("corrected.expt")[0]

    # Check that dx/dy lookup files have been set
    assert expt2.imageset.external_lookup.dx.filename
    assert expt2.imageset.external_lookup.dy.filename

    # Check mm positions are unaffected, but px positions are corrected
    mm1 = expt1.detector[0].get_beam_centre(expt1.beam.get_s0())
    mm2 = expt2.detector[0].get_beam_centre(expt2.beam.get_s0())
    px1 = expt1.detector[0].get_beam_centre_px(expt1.beam.get_s0())
    px2 = expt2.detector[0].get_beam_centre_px(expt2.beam.get_s0())
    assert mm1 == mm2
    assert px1[0] == px2[0] - 1
    assert px1[1] == px2[1] - 2


def test_ellipse_transforms():
    """See https://www.le.ac.uk/users/dsgp1/COURSES/TOPICS/quadrat.pdf for definitions"""

    # Generate random ellipse parameters
    phi = random.uniform(0, 360)
    l1 = random.uniform(0.5, 1.5)
    l2 = random.uniform(0.5, 1.5)

    # Check the transforms are inverses
    m1 = ellipse_to_circle_transform(phi, l1, l2)
    m2 = circle_to_ellipse_transform(phi, l1, l2)
    assert m2.elems == pytest.approx(m1.inverse().elems)

    # Generate some points around a circle
    ticks = np.arange(0, 2 * math.pi, math.pi / 10)
    radius = 5
    x = radius * np.cos(ticks)
    y = radius * np.sin(ticks)
    p1 = flex.vec2_double(zip(x, y))

    # Transform the points to an ellipse
    p2 = p1.__rmul__(m2)

    # Form the coefficients of the general ellipse equation,
    # a11 x^2 + 2 a12 xy + a22 y^2 = r^2
    cphi = math.cos(math.radians(phi))
    sphi = math.sin(math.radians(phi))
    a11 = l1 * cphi**2 + l2 * sphi**2
    a12 = (l2 - l1) * sphi * cphi
    a21 = a12
    a22 = l1 * sphi**2 + l2 * cphi**2
    A = matrix.sqr((a11, a12, a21, a22))

    # Calculate r^2 for the points on the ellipse, using the quadratic form x^T A x = r^2
    r_sq = p2.dot(p2.__rmul__(A))

    # Check that the points give the expected r^2
    assert r_sq == pytest.approx(radius**2)

    # When one of the scale parameters is 1.0, then the offsets should all have the same absolute
    # angle (the rotation fix from https://github.com/dials/dials/issues/3124#issuecomment-4109235695)
    l1 = 1.0
    m2 = circle_to_ellipse_transform(phi, l1, l2)
    p2 = p1.__rmul__(m2)
    offsets = p2 - p1
    x = flex.vec2_double(len(offsets), (1, 0))
    angles = abs(offsets.each_normalize().dot(x))
    assert angles.all_approx_equal(angles[0])


def test_elliptical_distortion_simple(run_in_tmp_path):
    """Create distortion maps for elliptical distortion using a dummy experiments
    with a small detector, for speed. Check those maps seem sensible"""

    # Make a detector model
    d = make_detector()

    # The beam is also essential for a experiments to be serialisable
    b = Beam((0, 0, 1), 1.0)

    # Create and write out a experiments
    imageset = ImageSet(ImageSetData(Reader(None, ["non-existent.cbf"]), None))
    imageset.set_detector(d)
    imageset.set_beam(b)
    experiments = ExperimentListFactory.from_imageset_and_crystal(imageset, None)
    experiments.as_json("dummy.expt")

    # Centre of distortion will be the far corner from the origin of the first
    # panel
    centre_xy = d[0].get_image_size_mm()

    # Generate distortion maps
    cmd = (
        "dials.generate_distortion_maps dummy.expt "
        "mode=ellipse centre_xy={},{} "
        "phi=90 l1=1.0 l2=0.95"
    ).format(*centre_xy)
    easy_run.fully_buffered(command=cmd).raise_if_errors()

    # Load the maps
    with open("dx.pickle", "rb") as f:
        dx = pickle.load(f)
    with open("dy.pickle", "rb") as f:
        dy = pickle.load(f)

    # Check there are 4 maps each
    assert len(dx) == len(dy) == 4

    # Ellipse has phi=90, so all correction is in the dx map. This is simplest
    # to understand because the fast direction is along the X axis.
    for arr in dy:
        assert min(arr) == pytest.approx(0.0)
        assert max(arr) == pytest.approx(0.0)

    # The ellipse correction is centred at the middle of the detector and all in
    # the X direction. Therefore we expect a few things from the dx maps:
    #
    # (1) Within each panel the columns of the array are identical.
    # (2) The two leftmost panels should be the same
    # (3) The two rightmost panels should be the same.
    # (4) One row from a panel on the left is a negated, reversed column from a
    #     panel on the right.
    #
    # All together expect the 4 dx maps to look something like this (indicative
    # values only):
    #
    # /-----------\ /-----------\
    # | 3  2  1  0| |-0 -1 -2 -3|
    # | 3  2  1  0| |-0 -1 -2 -3|
    # | 3  2  1  0| |-0 -1 -2 -3|
    # | 3  2  1  0| |-0 -1 -2 -3|
    # | 3  2  1  0| |-0 -1 -2 -3|
    # \-----------/ \-----------/
    # /-----------\ /-----------\
    # | 3  2  1  0| |-0 -1 -2 -3|
    # | 3  2  1  0| |-0 -1 -2 -3|
    # | 3  2  1  0| |-0 -1 -2 -3|
    # | 3  2  1  0| |-0 -1 -2 -3|
    # | 3  2  1  0| |-0 -1 -2 -3|
    # \-----------/ \-----------/

    # So the fundamental data is all in the first row of first panel's map
    dx0t = dx[0].matrix_transpose()
    row0 = dx0t.matrix_copy_column(0)

    # Test (1) from the above list for panel 0
    for i in range(1, d[0].get_image_size()[0]):
        assert row0 == pytest.approx(dx0t.matrix_copy_column(i))

    # Test (2)
    assert dx[0] == pytest.approx(dx[2])

    # Test (3)
    assert dx[1] == pytest.approx(dx[3])

    # Test (4)
    assert row0 == pytest.approx(
        -1.0 * dx[1].matrix_transpose().matrix_copy_column(0).reversed()
    )

    # Test (1) for panel 1 as well
    dx1t = dx[1].matrix_transpose()
    row0 = dx1t.matrix_copy_column(0)
    for i in range(1, d[0].get_image_size()[0]):
        assert row0 == pytest.approx(dx1t.matrix_copy_column(i))

    # The distortion is an ellipse pinched in by 5% (l2 = 0.95) horizontally.
    # Check that the distorted distance is 95% of the canonical distance for
    # the upper left pixel of the first panel.
    vec_centre_to_first_px = matrix.col(
        d[0].get_pixel_lab_coord((0.5, 0.5))
    ) - matrix.col(d[0].get_lab_coord(centre_xy))
    vec_distortion = matrix.col(
        (
            dx[0][0, 0] * d[0].get_pixel_size()[0],
            dy[0][0, 0] * d[0].get_pixel_size()[1],
            0.0,
        )
    )
    fast = matrix.col(d[0].get_fast_axis())
    assert (vec_centre_to_first_px + vec_distortion).dot(
        fast
    ) == 0.95 * vec_centre_to_first_px.dot(fast)


def create_distorted_ellipse_image(image_path, tmp_path, beam_centre_px, phi, l2):
    """Helper function to write an image consisting of 100 points around an ellipse"""

    subprocess.run(
        [
            shutil.which("dials.import"),
            image_path,
            f"fast_slow_beam_centre={beam_centre_px[0]},{beam_centre_px[1]}",
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    experiments = ExperimentList.from_file(tmp_path / "imported.expt")
    beam = experiments[0].beam
    panel = experiments[0].detector[0]

    # Put centre of distortion at the beam centre
    centre_xy = panel.get_beam_centre(beam.get_s0())

    # Get beam vector and two orthogonal vectors
    beamvec = matrix.col(beam.get_s0())
    bor1 = beamvec.ortho()
    bor2 = beamvec.cross(bor1)

    # Generate rays around a 2θ circle out to halfway to the panel edge
    d_min = panel.get_max_resolution_ellipse(beam)
    theta = math.asin(beam.get_wavelength() / (2 * d_min)) / 2
    n_rays = 100
    cone_base_centre = beamvec * math.cos(2.0 * theta)
    cone_base_radius = (beamvec * math.sin(2.0 * theta)).length()
    rad1 = bor1.normalize() * cone_base_radius
    rad2 = bor2.normalize() * cone_base_radius
    ticks = (2.0 * math.pi / n_rays) * flex.double_range(n_rays)
    offset1 = flex.vec3_double(n_rays, rad1) * flex.cos(ticks)
    offset2 = flex.vec3_double(n_rays, rad2) * flex.sin(ticks)
    rays = flex.vec3_double(n_rays, cone_base_centre) + offset1 + offset2

    # Get undistorted mm intersections on the detector
    circle_mm = flex.vec2_double((panel.get_ray_intersection(ray) for ray in rays))
    # from matplotlib import pyplot as plt
    # plt.scatter(*zip(*circle_mm))

    # Get the matrix to distort to the points to an ellipse
    l1 = 1.0
    m2 = circle_to_ellipse_transform(phi, l1, l2)

    # Distort the intersection points
    ellipse_mm = (circle_mm - centre_xy).__rmul__(m2) + centre_xy
    # plt.scatter(*zip(*ellipse_mm))
    # plt.gca().set_aspect("equal")
    # plt.show()

    # Get rays for the distorted intersections and get new intersections
    lab_coords = panel.get_lab_coord(ellipse_mm)
    rays = lab_coords.each_normalize() * (1.0 / beam.get_wavelength())
    panel = experiments[0].detector[0]
    intersections_px = flex.vec2_double(
        (panel.get_ray_intersection_px(ray) for ray in rays)
    )

    # Now write out a new image consisting only of the points around the ellipse
    with mrcfile.open(image_path, permissive=True) as mrc_in:
        data = mrc_in.data.copy()
        exttyp = mrc_in.header["exttyp"]
        extended_header = mrc_in.extended_header.copy()

    # Set the data to zero, then set the points around the ellipse to 100
    data[:] = 0
    for f, s in intersections_px:
        data[int(s), int(f)] = 100

    with mrcfile.new(tmp_path / "ellipse_001.mrc", overwrite=True) as mrc_out:
        mrc_out.set_data(data.astype(np.float32))
        mrc_out.header["exttyp"] = exttyp
        mrc_out.set_extended_header(extended_header)


def test_undistort_an_ellipse(dials_data, tmp_path):
    """Check that impact points around an ellipse in lab space on a simple
    detector can be undistorted into a circle in pixel space"""

    # Create an image containing spots in an ellipse with the minor axis 90% the
    # langth of the major axis, with the orientation of the major axis rotated by
    # a random angle from the fast axis. Shift the beam centre away from the image
    # centre for a more general test.
    beam_centre_px = (1040, 1010)
    phi = np.random.uniform(low=0.0, high=180.0)
    print(f"ellipse angle: φ={phi:.1f}°")
    l2 = 0.9
    image_path = dials_data("aluminium_standard") / "0p67_5s_0000.mrc"
    create_distorted_ellipse_image(image_path, tmp_path, beam_centre_px, phi, l2)

    # Import the new image and find the spots
    result = subprocess.run(
        [
            shutil.which("dials.import"),
            tmp_path / "ellipse_001.mrc",
            "fast_slow_beam_centre=1040,1010",
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    result = subprocess.run(
        [
            shutil.which("dials.find_spots"),
            tmp_path / "imported.expt",
            "min_spot_size=1",
        ],
        cwd=tmp_path,
        capture_output=True,
    )

    # Check that the the correct number of spots are found
    assert result.returncode == 0
    assert b"Saved 100 reflections to strong.refl" in result.stdout

    # Load the spots and calculate ellipse parameters as we would with the ellipse
    # tool in dials.image_viewer.
    strong = flex.reflection_table.from_file(tmp_path / "strong.refl")
    intersections_f, intersections_s, _ = strong["xyzobs.px.value"].parts()
    intersections_px = flex.vec2_double(intersections_f, intersections_s)
    try:
        ellipse = EllipseModel.from_estimate(np.array(intersections_px))
    except AttributeError:
        # Deprecated from skimage 0.26
        ellipse = EllipseModel()
        ellipse.estimate(np.array(intersections_px))
    phi_, a, b, centre_ = extract_ellipse_parameters(ellipse)
    l2_ = b / a
    assert phi_ == pytest.approx(phi, 0.1)
    assert l2_ == pytest.approx(l2, 0.01)
    assert centre_ == pytest.approx((beam_centre_px))

    # Generate and apply distortion maps using the ellipse parameters. Note the centre
    # must be in millimetres, not pixels.
    experiments = ExperimentList.from_file(tmp_path / "imported.expt")
    beam = experiments[0].beam
    panel = experiments[0].detector[0]
    centre_xy = panel.get_beam_centre(beam.get_s0())
    result = subprocess.run(
        [
            shutil.which("dials.generate_distortion_maps"),
            tmp_path / "imported.expt",
            "mode=ellipse",
            f"centre_xy={centre_xy[0]},{centre_xy[1]}",
            f"phi={phi_}",
            "l1=1.0",
            f"l2={l2_}",
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr

    result = subprocess.run(
        [
            shutil.which("dials.import"),
            tmp_path / "ellipse_001.mrc",
            "fast_slow_beam_centre=1040,1010",
            f"lookup.dx={tmp_path / 'dx.pickle'}",
            f"lookup.dy={tmp_path / 'dy.pickle'}",
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr

    # Now find spots again
    result = subprocess.run(
        [
            shutil.which("dials.find_spots"),
            tmp_path / "imported.expt",
            "min_spot_size=1",
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    assert result.returncode == 0
    strong = flex.reflection_table.from_file(tmp_path / "strong.refl")

    # The pixel positions should be unchanged.
    intersections_f_, intersections_s_, _ = strong["xyzobs.px.value"].parts()
    assert intersections_f_.all_eq(intersections_f)
    assert intersections_s_.all_eq(intersections_s)
    intersections_px = flex.vec2_double(intersections_f_, intersections_s_)

    # Now calculate mm positions, which should be undistorted into a circle
    experiments = ExperimentList.from_file(tmp_path / "imported.expt")
    beam = experiments[0].beam
    panel = experiments[0].detector[0]
    intersections_mm = panel.pixel_to_millimeter(intersections_px)

    # Get the radius of each of the intersections and check fractional error.
    # With high distortion of l2=0.9 we found error as high as 1.2%
    shifted = intersections_mm - centre_xy
    x, y = shifted.parts()
    radius = flex.sqrt(x * x + y * y)
    mm_error = (max(radius) - min(radius)) / flex.mean(radius)
    print(f"mm_error = {mm_error * 100:.1f}%")
    assert mm_error < 0.013

    # Wiht l2=0.9 we seem to get radial errors of up to around 3 pixels. With
    # l2=0.95 the error is up to about 1.1 pixels. This might be acceptable.
    assert radius.as_numpy_array() == pytest.approx(
        flex.mean(radius), abs=3.1 * panel.get_pixel_size()[0]
    )

    # Visually inspect the mm corrected positions
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()
    intersections_x, intersections_y = intersections_mm.parts()
    ax.scatter(intersections_x, intersections_y)
    ax.yaxis.set_inverted(True)
    ax.set_xlim(0, 2048 * panel.get_pixel_size()[0])
    ax.set_ylim(0, 2048 * panel.get_pixel_size()[1])
    plt.gca().set_aspect("equal", adjustable="box")
    plt.title("mm corrected positions")
    plt.show()

    # Load the experiment and extract the correction maps
    experiments = ExperimentList.from_file(tmp_path / "imported.expt")
    dx = experiments[0].imageset.external_lookup.dx.data.tile(0).data()
    dy = experiments[0].imageset.external_lookup.dy.data.tile(0).data()

    # Check that the correction at the beam centre is very small. It isn't
    # actually the smallest magnitude correction across the whole image. There
    # is a line of pixels with close to zero correction from one side of the
    # image to the other, passing through the beam centre. We'll accept that
    # the correction at the pixel containing the beam centre should be less than
    # 1/5 of a pixel.
    assert abs(dx[int(beam_centre_px[1]), int(beam_centre_px[0])]) < 0.2
    assert abs(dy[int(beam_centre_px[1]), int(beam_centre_px[0])]) < 0.2

    # For each intersection get the correction encoded by the distortion maps.
    # In the pixel-to-millimetre transform, the dx and dy values are first
    # *subtracted* from the pixel positions and then the result is converted
    # to millimetre via the parallax correction.
    corrections_f = []
    corrections_s = []
    for spot in intersections_px:
        index = (int(spot[1]), int(spot[0]))
        corrections_f.append(dx[*index])
        corrections_s.append(dy[*index])
    corrections_f = flex.double(corrections_f)
    corrections_s = flex.double(corrections_s)

    corrected_intersections_f = intersections_f - corrections_f
    corrected_intersections_s = intersections_s - corrections_s

    from matplotlib import pyplot as plt

    fig, ax = plt.subplots()
    ax.scatter(intersections_f, intersections_s)
    ax.scatter(corrected_intersections_f, corrected_intersections_s)
    ax.quiver(
        corrected_intersections_f,
        corrected_intersections_s,
        corrections_f,
        corrections_s,
        angles="xy",
        scale_units="xy",
        scale=1,
        color="red",
    )
    ax.yaxis.set_inverted(True)
    plt.gca().set_aspect("equal", adjustable="box")
    plt.title("pixel intersections before/after correction")
    plt.show()

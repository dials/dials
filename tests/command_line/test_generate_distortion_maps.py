from __future__ import annotations

import math
import os
import pickle
import random
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest

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


def test_translate(dials_regression: Path, run_in_tmp_path):
    """Test as written in https://github.com/dials/dials/issues/471. This
    is pretty slow!"""

    # use the i04_weak_data for this test
    data_dir = os.path.join(dials_regression, "image_examples", "DLS_I04")
    image_path = os.path.join(data_dir, "grid_full_cbf_0005.cbf")

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
        "phi=0 l1=1.0 l2=0.95"
    ).format(*centre_xy)
    easy_run.fully_buffered(command=cmd).raise_if_errors()

    # Load the maps
    with open("dx.pickle", "rb") as f:
        dx = pickle.load(f)
    with open("dy.pickle", "rb") as f:
        dy = pickle.load(f)

    # Check there are 4 maps each
    assert len(dx) == len(dy) == 4

    # Ellipse has phi=0, so all correction is in the dy map
    for arr in dx:
        assert min(arr) == max(arr) == 0.0

    # The ellipse correction is centred at the middle of the detector and all in
    # the Y direction. Therefore we expect a few things from the dy maps:
    #
    # (1) Within each panel the columns of the array are identical.
    # (2) The two upper panels should be the same
    # (3) The two lower panels should be the same.
    # (4) One column from an upper panel is a negated, reversed column from a
    #     lower panel.
    #
    # All together expect the 4 dy maps to look something like this:
    #
    # /-----------\ /-----------\
    # |-4 -4 -4 -4| |-4 -4 -4 -4|
    # |-3 -3 -3 -3| |-3 -3 -3 -3|
    # |-2 -2 -2 -2| |-2 -2 -2 -2|
    # |-1 -1 -1 -1| |-1 -1 -1 -1|
    # | 0  0  0  0| | 0  0  0  0|
    # \-----------/ \-----------/
    # /-----------\ /-----------\
    # | 0  0  0  0| | 0  0  0  0|
    # | 1  1  1  1| | 1  1  1  1|
    # | 2  2  2  2| | 2  2  2  2|
    # | 3  3  3  3| | 3  3  3  3|
    # | 4  4  4  4| | 4  4  4  4|
    # \-----------/ \-----------/

    # So the fundamental data is all in the first column of first panel's map
    col0 = dy[0].matrix_copy_column(0)

    # The correction should be 5% of the distance from the ellipse centre to a
    # corrected pixel (l2 = 0.95 above) along the slow axis. Check that is the
    # case (for the first pixel at least)
    vec_centre_to_first_px = matrix.col(
        d[0].get_pixel_lab_coord((0.5, 0.5))
    ) - matrix.col(d[0].get_lab_coord(centre_xy))
    dist_centre_to_first_px = vec_centre_to_first_px.dot(
        matrix.col(d[0].get_slow_axis())
    )
    corr_mm = dist_centre_to_first_px * 0.05
    corr_px = corr_mm / d[0].get_pixel_size()[1]
    assert col0[0] == pytest.approx(corr_px)

    # Test (1) from above list for panel 0
    for i in range(1, d[0].get_image_size()[0]):
        assert (col0 == dy[0].matrix_copy_column(i)).all_eq(True)

    # Test (2)
    assert (dy[0] == dy[1]).all_eq(True)

    # Test (3)
    assert (dy[2] == dy[3]).all_eq(True)

    # Test (4)
    assert col0 == pytest.approx(-1.0 * dy[2].matrix_copy_column(0).reversed())

    # Test (1) for panel 2 as well, which then covers everything needed
    col0 = dy[2].matrix_copy_column(0)
    for i in range(1, d[0].get_image_size()[0]):
        assert (col0 == dy[2].matrix_copy_column(i)).all_eq(True)


def test_undistort_an_ellipse(dials_data, tmp_path):
    """Check that impact points around an ellipse in lab space on a simple
    detector can be undistorted into a circle in pixel space"""

    # Use a single-panel 3D ED image for this
    image_path = (
        dials_data("image_examples", pathlib=True) / "TIMEPIX_SU_516-stdgoni_0001.img"
    )
    result = subprocess.run(
        [
            shutil.which("dials.import"),
            image_path,
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    experiments = ExperimentList.from_file(tmp_path / "imported.expt")
    beam = experiments[0].beam
    panel = experiments[0].detector[0]

    # Put centre of distortion at the beam centre
    centre_xy = panel.get_beam_centre(beam.get_s0())
    centre_px = panel.millimeter_to_pixel(centre_xy)

    # Get beam vector and two orthogonal vectors
    beamvec = matrix.col(beam.get_s0())
    bor1 = beamvec.ortho()
    bor2 = beamvec.cross(bor1)

    # Generate rays at a 2θ circle out to halfway to the panel edge
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

    # Get the matrix to distort to a rotated ellipse
    phi = 15
    l1 = 1.0
    l2 = 0.95
    m2 = circle_to_ellipse_transform(phi, l1, l2)

    # Distort the intersection points
    ellipse_mm = (circle_mm - centre_xy).__rmul__(m2) + centre_xy
    # plt.scatter(*zip(*ellipse_mm))
    # plt.gca().set_aspect("equal")
    # plt.show()

    # Get rays for the distorted intersections
    lab_coords = panel.get_lab_coord(ellipse_mm)
    rays = lab_coords.each_normalize() * (1.0 / beam.get_wavelength())

    # Generate and apply distortion maps
    result = subprocess.run(
        [
            shutil.which("dials.generate_distortion_maps"),
            tmp_path / "imported.expt",
            "mode=ellipse",
            f"centre_xy={centre_xy[0]},{centre_xy[1]}",
            f"phi={phi}",
            f"l1={l1}",
            f"l2={l2}",
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr

    result = subprocess.run(
        [
            shutil.which("dials.import"),
            image_path,
            f"lookup.dx={tmp_path / 'dx.pickle'}",
            f"lookup.dy={tmp_path / 'dy.pickle'}",
        ],
        cwd=tmp_path,
        capture_output=True,
    )
    assert not result.returncode and not result.stderr

    # Load the experiment with correction maps and calculate ray intersections
    experiments = ExperimentList.from_file(tmp_path / "imported.expt")
    panel = experiments[0].detector[0]
    intersections_px = flex.vec2_double(
        (panel.get_ray_intersection_px(ray) for ray in rays)
    )

    # Check that the pixel intersections really are circular
    shifted = intersections_px - centre_px
    x, y = shifted.parts()
    r = flex.sqrt(x * x + y * y)

    # Seem to have errors of half a pixel or so...
    assert r.as_numpy_array() == pytest.approx(flex.mean(r), abs=0.5)

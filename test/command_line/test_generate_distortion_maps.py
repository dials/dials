from __future__ import absolute_import, division, print_function

import pytest
from dxtbx.model.detector import Detector
from dxtbx.model.beam import Beam
from scitbx import matrix
from dxtbx.imageset import ImageSet
from dxtbx.imageset import ImageSetData
from dxtbx.format.Format import Reader
from dxtbx.model.experiment_list import ExperimentListFactory

from libtbx import easy_run
import six.moves.cPickle as pickle
import os


def make_detector():
    """Make a dummy 4 panel detector with not many pixels to ensure test runs
    quickly"""
    pixel_size_x = 0.1
    pixel_size_y = 0.1
    npixels_per_panel_x = 50
    npixels_per_panel_y = 50
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


@pytest.mark.xfail(reason="dx,dy maps not loaded from json")
def test_translate(dials_regression, run_in_tmpdir):
    """Test as written in https://github.com/dials/dials/issues/471. This
    is pretty slow!"""

    # use the i04_weak_data for this test
    data_dir = os.path.join(dials_regression, "image_examples", "DLS_I04")
    image_path = os.path.join(data_dir, "grid_full_cbf_0005.cbf")

    # Generate distortion maps
    cmd = ("dials.generate_distortion_maps {0} " "dx=1 dy=2").format(image_path)
    easy_run.fully_buffered(command=cmd).raise_if_errors()

    # Import without correction
    cmd = ("dials.import {0}").format(image_path)
    easy_run.fully_buffered(command=cmd).raise_if_errors()
    # experiments = ExperimentListFactory.from_serialized_format("imported.expt")[0]
    # det1 = experiments.detectors()[0]

    # Import with correction
    cmd = (
        "dials.import {0} dx=dx.pickle dy=dy.pickle "
        "output.experiments=corrected.expt"
    ).format(image_path)
    easy_run.fully_buffered(command=cmd).raise_if_errors()
    # experiments2 = ExperimentListFactory.from_serialized_format("corrected.expt")[0]
    # det2 = experiments2.detectors()[0]

    # FIXME, why doesn't db2 have dx, dy set?
    assert db2.extract_imagesets()[0].external_lookup.dx.filename

    # FIXME finish test by comparing px to mm positions for det1, det2


def test_elliptical_distortion(run_in_tmpdir):
    """Create distortion maps for elliptical distortion using a dummy experiments
    with a small detector, for speed. Check those maps seem sensible"""

    # Make a detector model
    d = make_detector()

    # The beam is also essential for a experiments to be serialisable
    b = Beam((0, 0, 1), 1.0)

    # Create and write out a experiments
    imageset = ImageSet(ImageSetData(Reader(["non-existent.cbf"]), None))
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
    for i in range(1, 50):
        assert (col0 == dy[0].matrix_copy_column(i)).all_eq(True)

    # Test (2)
    assert (dy[0] == dy[1]).all_eq(True)

    # Test (3)
    assert (dy[2] == dy[3]).all_eq(True)

    # Test (4)
    assert col0 == pytest.approx(-1.0 * dy[2].matrix_copy_column(0).reversed())

    # Test (1) for panel 2 as well, which then covers everything needed
    col0 = dy[2].matrix_copy_column(0)
    for i in range(1, 50):
        assert (col0 == dy[2].matrix_copy_column(i)).all_eq(True)

from __future__ import absolute_import, division, print_function

import pytest
from dials.command_line.frame_orientations import extract_experiment_data


def test_extract_experiment_data():
    """Test basic operation of the extract_experiment_data function. Does not
    test extraction of data from scan-varying models"""

    # Set up an Experiment with idealised geometry
    from dxtbx.model import BeamFactory
    from dxtbx.model import GoniometerFactory
    from dxtbx.model import Crystal
    from dxtbx.model import ScanFactory
    from dxtbx.model.experiment_list import Experiment

    beam = BeamFactory.make_beam(unit_s0=(0, 0, -1), wavelength=1.0)
    goniometer = GoniometerFactory.known_axis((1, 0, 0))
    a = (100, 0, 0)
    b = (0, 90, 0)
    c = (0, 0, 80)
    crystal = Crystal(a, b, c, space_group_symbol="P1")
    scan = ScanFactory.make_scan(
        image_range=(1, 91),
        exposure_times=0.1,
        oscillation=(0, 1.0),
        epochs=list(range(91)),
        deg=True,
    )

    exp = Experiment(beam=beam, goniometer=goniometer, scan=scan, crystal=crystal)

    # Extract experiment data
    dat = extract_experiment_data(exp, scale=100)

    # Check results are as expected
    za = dat["zone_axes"]

    # At the first image the c axis is aligned antiparallel with the beam vector,
    # while the a and b axes are orthogonal. The zone axis calculation is scaled
    # by 100 (i.e. the max cell dimension, which is the default), therefore we
    # expect the zone axis [uvw] = [0 0 -100/80]
    assert za[0].elems == pytest.approx((0, 0, -100 / 80))

    # At the start of the 91st image the crystal has rotated by 90 degrees, so
    # now c is orthogonal to the beam while b is anti-parallel to it. The zone
    # axis is now expected to be [uvw] = [0 -100/90 0]
    assert za[-1].elems == pytest.approx((0, -100 / 90, 0))

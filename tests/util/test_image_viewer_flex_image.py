"""The image viewer composites multi-panel detectors onto one picture. Panels
with different pixel sizes and different pixel counts must land at a
consistent scale and position, and picture pixels must look up the right
panel."""

from __future__ import annotations

import pytest

from dxtbx.model import Beam, Detector, Panel
from scitbx.array_family import flex

from dials.util.image_viewer.slip_viewer.flex_image import (
    get_flex_image_multipanel,
    get_picture_pixel_size,
)


def _detector(specs):
    """specs: list of (pixel_size_mm, (nfast, nslow), origin_mm). All panels
    lie in the plane z = -100 with fast = +x, slow = -y."""
    det = Detector()
    data = []
    for i, (px, (nf, ns), origin) in enumerate(specs):
        p = Panel()
        p.set_name("p%d" % i)
        p.set_image_size((nf, ns))
        p.set_pixel_size((px, px))
        p.set_frame((1, 0, 0), (0, -1, 0), origin)
        p.set_trusted_range((-1, 1e6))
        det.add_panel(p)
        arr = flex.double(flex.grid(ns, nf), float(i + 1))
        data.append(arr)
    return det, data


@pytest.mark.parametrize(
    "specs",
    [
        # 2 x 2 quadrants, 20 x 20 mm each with 1 mm gaps, mixed pixel sizes
        [
            (0.1, (200, 200), (-20.5, 20.5, -100)),
            (0.05, (400, 400), (0.5, 20.5, -100)),
            (0.2, (100, 100), (-20.5, -0.5, -100)),
            (0.1, (200, 200), (0.5, -0.5, -100)),
        ],
        # same pixel size, unequal panel sizes
        [
            (0.1, (300, 250), (-20.5, 20.5, -100)),
            (0.1, (100, 250), (10.5, 20.5, -100)),
            (0.1, (300, 150), (-20.5, -5.5, -100)),
            (0.1, (100, 150), (10.5, -5.5, -100)),
        ],
    ],
)
def test_multipanel_picture_is_consistent(specs):
    det, data = _detector(specs)
    beam = Beam()
    beam.set_unit_s0((0, 0, -1))
    beam.set_wavelength(1.0)

    fi = get_flex_image_multipanel(det, data, beam)
    picture_px = get_picture_pixel_size(det)
    assert picture_px == min(s[0] for s in specs)

    # The same lab position must map to the same picture position whichever
    # panel it is taken from: compare panel corners that touch across a gap.
    def picture_of_lab(k, f, s):
        return fi.tile_readout_to_picture(k, s, f)

    # extents in picture pixels scale with pixel size
    for k, panel in enumerate(det):
        nf, ns = panel.get_image_size()
        a = picture_of_lab(k, 0, 0)
        b = picture_of_lab(k, nf, ns)
        scale = panel.get_pixel_size()[0] / picture_px
        assert abs(abs(b[0] - a[0]) - ns * scale) < 1e-6
        assert abs(abs(b[1] - a[1]) - nf * scale) < 1e-6

    # lab -> picture is one affine map for all panels: fit on panel 0, test all
    import numpy as np

    def samples(k):
        p = det[k]
        nf, ns = p.get_image_size()
        rows = []
        for s in (0, ns // 2, ns):
            for f in (0, nf // 2, nf):
                pic = picture_of_lab(k, f, s)
                lab = p.get_pixel_lab_coord((f, s))
                rows.append((pic[0], pic[1], lab[0], lab[1]))
        return np.array(rows)

    S0 = samples(0)
    A, *_ = np.linalg.lstsq(
        np.c_[S0[:, 2], S0[:, 3], np.ones(len(S0))], S0[:, :2], rcond=None
    )
    for k in range(len(det)):
        S = samples(k)
        pred = np.c_[S[:, 2], S[:, 3], np.ones(len(S))] @ A
        assert np.abs(pred - S[:, :2]).max() < 1e-6

    # picture -> readout lookups return the panel the point lies on
    fi.setWindow(0, 0, 1)
    fi.adjust(color_scheme=0)
    fi.prep_string()
    for k, panel in enumerate(det):
        nf, ns = panel.get_image_size()
        pc = picture_of_lab(k, nf // 2, ns // 2)
        z = fi.picture_to_readout(pc[0], pc[1])
        assert int(z[2]) == k
        assert abs(z[0] - ns // 2) < 1 and abs(z[1] - nf // 2) < 1
        # a point just outside the panel's far corner is in a gap or outside
        corner = picture_of_lab(k, nf, ns)
        direction = np.sign(np.array(corner) - np.array(pc))
        outside = np.array(corner) + 3 * direction
        z_out = fi.picture_to_readout(outside[0], outside[1])
        assert int(z_out[2]) != k

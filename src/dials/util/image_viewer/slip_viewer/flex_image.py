from __future__ import annotations

import math

import scitbx.matrix
from dxtbx.model.detector_helpers import get_panel_projection_2d_from_axes
from scitbx.array_family import flex


def get_flex_image(
    data,
    vendortype,
    binning=1,
    brightness=1.0,
    saturation=65535.0,
    show_untrusted=False,
    color_scheme=0,
):
    # This is a combination of the get_data_type() and get_flex_image()
    # functions from iotbx.detectors.detectorbase.  XXX This may turn
    # out to be generally useful (see
    # e.g. rstbx.viewer.create_flex_image()), but where to place it?
    # dxtbx Format class?
    typehash = str(data.__class__)
    if typehash.find("int") >= 0:
        from iotbx.detectors import FlexImage
    elif typehash.find("double") >= 0:
        from iotbx.detectors import FlexImage_d as FlexImage

    return FlexImage(
        binning=binning,
        brightness=brightness,
        rawdata=data,
        saturation=int(round(saturation)),
        vendortype=vendortype,
        show_untrusted=show_untrusted,
        color_scheme=color_scheme,
    )


def get_flex_image_multipanel(
    detector,
    image_data,
    beam,
    brightness=1.0,
    binning=1,
    show_untrusted=False,
    color_scheme=0,
):
    # From xfel.cftbx.cspad_detector.readHeader() and
    # xfel.cftbx.cspad_detector.get_flex_image().  XXX Is it possible to
    # merge this with get_flex_image() above?  XXX Move to dxtbx Format
    # class (or a superclass for multipanel images)?

    from iotbx.detectors import generic_flex_image
    from libtbx.test_utils import approx_equal

    assert len(detector) == len(image_data), (len(detector), len(image_data))

    # Determine next multiple of eight of the largest panel size.
    data_max_focus = None
    for data in image_data:
        if data_max_focus is None:
            data_max_focus = data.focus()
        else:
            data_max_focus = (
                max(data_max_focus[0], data.focus()[0]),
                max(data_max_focus[1], data.focus()[1]),
            )
    data_padded = (
        8 * int(math.ceil(data_max_focus[0] / 8)),
        8 * int(math.ceil(data_max_focus[1] / 8)),
    )

    # Assert that all saturated values are equal and not None.  While
    # dxtbx records a separated trusted_range for each panel,
    # generic_flex_image supports only accepts a single common value for
    # the saturation.
    saturation = None
    for panel in detector:
        if saturation is None:
            saturation = panel.get_trusted_range()[1]
        else:
            assert approx_equal(saturation, panel.get_trusted_range()[1])
    assert saturation is not None

    # Create rawdata and flex_image_multipanel before populating it.
    rawdata = flex.double(flex.grid(len(detector) * data_padded[0], data_padded[1]))
    flex_image_multipanel = generic_flex_image(
        rawdata=rawdata,
        binning=binning,
        size1_readout=data_max_focus[0],
        size2_readout=data_max_focus[1],
        brightness=brightness,
        saturation=saturation,
        show_untrusted=show_untrusted,
        color_scheme=color_scheme,
    )

    # Calculate the average beam center across all panels, in meters
    # not sure this makes sense for detector which is not on a plane?
    beam_center = scitbx.matrix.col((0, 0, 0))
    npanels = 0
    for panel in detector:
        try:
            beam_center += scitbx.matrix.col(panel.get_beam_centre_lab(beam.get_s0()))
            npanels += 1
        except RuntimeError:  # catch DXTBX_ASSERT for no intersection
            pass
    beam_center /= npanels / 1e-3

    # XXX If a point is contained in two panels simultaneously, it will
    # be assigned to the panel defined first.  XXX Use a Z-buffer
    # instead?
    for i, panel in enumerate(detector):

        # Determine the pixel size for the panel (in meters), as pixel
        # sizes need not be identical.
        data = image_data[i]

        rawdata.matrix_paste_block_in_place(
            block=data.as_double(), i_row=i * data_padded[0], i_column=0
        )

        # If the panel already has a 2d projection then use it
        if panel.get_projection_2d():
            panel_r, panel_t = panel.get_projection_2d()
        else:
            if getattr(detector, "projection", "lab") == "image":
                # Get axes from precalculated 2D projection.
                origin_2d, fast_2d, slow_2d = detector.projection_2d_axes
                fast = scitbx.matrix.col(fast_2d[i] + (0,))
                slow = scitbx.matrix.col(slow_2d[i] + (0,))
                origin = scitbx.matrix.col(origin_2d[i] + (0,)) * 1e-3
            else:
                # Get unit vectors in the fast and slow directions, as well as the
                # the locations of the origin and the center of the panel, in
                # meters. The origin is taken w.r.t. to average beam center of all
                # panels. This avoids excessive translations that can result from
                # rotations around the laboratory origin. Related to beam centre above
                # and dials#380 not sure this is right for detectors which are not
                # coplanar since system derived from first panel...
                fast = scitbx.matrix.col(panel.get_fast_axis())
                slow = scitbx.matrix.col(panel.get_slow_axis())
                origin = scitbx.matrix.col(panel.get_origin()) * 1e-3 - beam_center

            panel_r, panel_t = get_panel_projection_2d_from_axes(
                panel, data, fast, slow, origin
            )

        flex_image_multipanel.add_transformation_and_translation(panel_r, panel_t)

    flex_image_multipanel.followup_brightness_scale()
    return flex_image_multipanel

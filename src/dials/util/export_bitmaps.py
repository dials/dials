from __future__ import annotations

import enum
import logging
import math
import sys
from collections.abc import Iterator, Sequence

from PIL import Image, ImageDraw, ImageFont

from cctbx import miller, sgtbx, uctbx
from dxtbx.model import ImageSet
from dxtbx.model.detector_helpers import get_detector_projection_2d_axes
from iotbx.detectors import FlexImage, FlexImage_d

from dials.algorithms.image.threshold import DispersionThresholdDebug
from dials.array_family import flex
from dials.util.image_viewer import calculate_isoresolution_lines
from dials.util.image_viewer.slip_viewer.flex_image import (
    get_flex_image,
    get_flex_image_multipanel,
)

logger = logging.getLogger(__name__)

HEXAGONAL_ICE_UNIT_CELL = uctbx.unit_cell((4.498, 4.498, 7.338, 90, 90, 120))
HEXAGONAL_ICE_SPACE_GROUP = sgtbx.space_group_info(194).group()


class Projection(enum.Enum):
    LAB = "lab"
    IMAGE = "image"


class Display(enum.Enum):
    IMAGE = "image"
    MEAN = "mean"
    VARIANCE = "variance"
    DISPERSION = "dispersion"
    SIGMA_B = "sigma_b"
    SIGMA_S = "sigma_s"
    THRESHOLD = "threshold"
    GLOBAL_THRESHOLD = "global_threshold"


class ColourScheme(enum.Enum):
    GREYSCALE = 0
    RAINBOW = 1
    HEATMAP = 2
    INVERSE_GREYSCALE = 3


def imageset_as_flex_image(
    imageset: ImageSet,
    images: Sequence[int],
    brightness: float = 100,
    binning: int = 1,
    projection: Projection = Projection.IMAGE,
    saturation: int = 0,
    show_mask: bool = False,
    display: Display = Display.IMAGE,
    colour_scheme: ColourScheme = ColourScheme.GREYSCALE,
    gain: float = 1,
    nsigma_b: float = 6,
    nsigma_s: float = 3,
    global_threshold: float = 0,
    min_local: int = 2,
    kernel_size: tuple[int, int] = (3, 3),
) -> Iterator[FlexImage | FlexImage_d]:
    brightness = brightness / 100
    # check that binning is a power of 2
    if not (binning > 0 and ((binning & (binning - 1)) == 0)):
        raise ValueError("binning must be a power of 2")

    detector = imageset.get_detector()

    # Furnish detector with 2D projection axes
    detector.projected_2d = get_detector_projection_2d_axes(detector)
    detector.projection = projection.value

    # XXX is this inclusive or exclusive?
    saturation = detector[0].get_trusted_range()[1]
    scan = imageset.get_scan()
    if scan is not None and not scan.is_still():  # and not images:
        start, _ = scan.get_image_range()
    else:
        start = 1

    # If the user specified an image range index, only export those
    n_images = len(imageset)
    for i_image in images:
        if (i_image < start) or (i_image >= start + n_images):
            raise ValueError(
                f"Image {i_image} outside of scan range {start},{start + n_images - 1}"
            )
        image = imageset.get_raw_data(i_image - start)

        mask = imageset.get_mask(i_image - start)
        if mask is None:
            mask = [p.get_trusted_range_mask(im) for im, p in zip(image, detector)]

        if show_mask:
            for rd, m in zip(image, mask):
                rd.set_selected(~m, -2)

        image = image_filter(
            image,
            mask,
            display=display,
            gain_value=gain,
            nsigma_b=nsigma_b,
            nsigma_s=nsigma_s,
            global_threshold=global_threshold,
            min_local=min_local,
            kernel_size=kernel_size,
        )

        if len(detector) > 1:
            # FIXME This doesn't work properly, as flex_image.size2() is incorrect
            # also binning doesn't work
            flex_image = get_flex_image_multipanel(
                brightness=brightness,
                detector=detector,
                image_data=image,
                binning=binning,
                beam=imageset.get_beam(),
                show_untrusted=show_mask,
            )
        else:
            flex_image = get_flex_image(
                brightness=brightness,
                data=image[0],
                binning=binning,
                saturation=saturation,
                vendortype="nonsense",
                show_untrusted=show_mask,
            )

        flex_image.setWindow(0, 0, 1)
        flex_image.adjust(color_scheme=colour_scheme.value)

        # now export as a bitmap
        flex_image.prep_string()

        yield flex_image


def image_filter(
    raw_data,
    mask,
    display: Display,
    gain_value,
    nsigma_b,
    nsigma_s,
    global_threshold,
    min_local,
    kernel_size,
):
    if display is Display.IMAGE:
        return raw_data

    assert gain_value > 0
    gain_map = [flex.double(rd.accessor(), gain_value) for rd in raw_data]

    kabsch_debug_list = [
        DispersionThresholdDebug(
            data.as_double(),
            mask[i_panel],
            gain_map[i_panel],
            kernel_size,
            nsigma_b,
            nsigma_s,
            global_threshold,
            min_local,
        )
        for i_panel, data in enumerate(raw_data)
    ]

    if display is Display.MEAN:
        display_data = [kabsch.mean() for kabsch in kabsch_debug_list]
    elif display is Display.VARIANCE:
        display_data = [kabsch.variance() for kabsch in kabsch_debug_list]
    elif display is Display.DISPERSION:
        display_data = [kabsch.index_of_dispersion() for kabsch in kabsch_debug_list]
    elif display is Display.SIGMA_B:
        cv = [kabsch.index_of_dispersion() for kabsch in kabsch_debug_list]
        display_data = (kabsch.cv_mask() for kabsch in kabsch_debug_list)
        display_data = [_mask.as_1d().as_double() for _mask in display_data]
        for i, _mask in enumerate(display_data):
            _mask.reshape(cv[i].accessor())
    elif display is Display.SIGMA_S:
        cv = [kabsch.index_of_dispersion() for kabsch in kabsch_debug_list]
        display_data = (kabsch.value_mask() for kabsch in kabsch_debug_list)
        display_data = [_mask.as_1d().as_double() for _mask in display_data]
        for i, _mask in enumerate(display_data):
            _mask.reshape(cv[i].accessor())
    elif display is Display.GLOBAL_THRESHOLD:
        cv = [kabsch.index_of_dispersion() for kabsch in kabsch_debug_list]
        display_data = (kabsch.global_mask() for kabsch in kabsch_debug_list)
        display_data = [_mask.as_1d().as_double() for _mask in display_data]
        for i, _mask in enumerate(display_data):
            _mask.reshape(cv[i].accessor())
    elif display is Display.THRESHOLD:
        cv = [kabsch.index_of_dispersion() for kabsch in kabsch_debug_list]
        display_data = (kabsch.final_mask() for kabsch in kabsch_debug_list)
        display_data = [_mask.as_1d().as_double() for _mask in display_data]
        for i, _mask in enumerate(display_data):
            _mask.reshape(cv[i].accessor())

    return display_data


def _draw_resolution_rings_impl(
    imageset: ImageSet,
    pil_image: Image.Image,
    flex_image: FlexImage | FlexImage_d,
    spacings: flex.double,
    fill: str,
    fontsize: int | None,
    binning: int,
):
    segments, res_labels = calculate_isoresolution_lines(
        spacings,
        imageset.get_beam(),
        imageset.get_detector(),
        flex_image,
        binning=binning,
    )
    draw = ImageDraw.Draw(pil_image)
    for segment in segments:
        draw.line(segment, fill=fill, width=2)
    if fontsize:
        try:
            # Only import matplotlib if we absolutely need to, and don't
            # override the backend if the user has already imported
            if "matplotlib" not in sys.modules:
                import matplotlib

                matplotlib.use("Agg")
            from matplotlib.font_manager import FontProperties, fontManager

            font_filename = fontManager.findfont(
                FontProperties(size=fontsize / binning**0.5), fontext="ttf"
            )
            font = ImageFont.truetype(
                font_filename,
                size=math.ceil(fontsize / binning**0.5),
            )
        except (ImportError, ValueError):
            # Revert to default bitmap font if we must, but fontsize will not work
            logger.warning("Could not find default font, using fallback")
            font = ImageFont.load_default()

        for x, y, label in res_labels:
            draw.text((x, y), label, fill=fill, font=font)


def draw_resolution_rings(
    imageset: ImageSet,
    pil_image: Image.Image,
    flex_image: FlexImage | FlexImage_d,
    n_rings: int = 5,
    spacings: list | None = None,
    fill: str = "red",
    fontsize: int | None = 30,
    binning: int = 1,
):
    if not spacings:
        d_min = imageset.get_detector().get_max_resolution(imageset.get_beam().get_s0())
        d_star_sq_max = uctbx.d_as_d_star_sq(d_min)

        step = d_star_sq_max / (n_rings + 1)
        spacings = flex.double(
            [uctbx.d_star_sq_as_d((i + 1) * step) for i in range(0, n_rings)]
        )
    else:
        spacings = flex.double(spacings)

    _draw_resolution_rings_impl(
        imageset,
        pil_image,
        flex_image,
        spacings,
        fill,
        fontsize,
        binning,
    )


def draw_ice_rings(
    imageset: ImageSet,
    pil_image: Image.Image,
    flex_image: FlexImage | FlexImage_d,
    unit_cell: uctbx.unit_cell = HEXAGONAL_ICE_UNIT_CELL,
    space_group: sgtbx.space_group = HEXAGONAL_ICE_SPACE_GROUP,
    fill: str = "blue",
    fontsize: int | None = 30,
    binning: int = 1,
):
    d_min = imageset.get_detector().get_max_resolution(imageset.get_beam().get_s0())
    unit_cell = space_group.average_unit_cell(unit_cell)
    generator = miller.index_generator(unit_cell, space_group.type(), False, d_min)
    indices = generator.to_array()
    spacings = flex.sorted(unit_cell.d(indices))
    _draw_resolution_rings_impl(
        imageset,
        pil_image,
        flex_image,
        spacings,
        fill,
        fontsize,
        binning,
    )

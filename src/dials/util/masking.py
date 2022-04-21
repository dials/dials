from __future__ import annotations

import copy
import logging
import math
import time
from collections import namedtuple
from typing import Tuple

import libtbx.phil
from cctbx import crystal
from dxtbx.masking import (
    mask_untrusted_circle,
    mask_untrusted_polygon,
    mask_untrusted_rectangle,
)
from dxtbx.model import ImageSet, SimplePxMmStrategy
from iotbx.phil import parse

from dials.array_family import flex
from dials.util.ext import ResolutionMaskGenerator

logger = logging.getLogger(__name__)

phil_scope = parse(
    """
  border = 0
    .type = int
    .help = "The border around the edge of the image."

  d_min = None
    .help = "The high resolution limit in Angstrom for a pixel to be"
            "accepted by the filtering algorithm."
    .type = float(value_min=0)

  d_max = None
    .help = "The low resolution limit in Angstrom for a pixel to be"
            "accepted by the filtering algorithm."
    .type = float(value_min=0)

  disable_parallax_correction = False
    .help = "Set to ``True`` to use a faster, but less accurate, simple px-to-mm "
            "mapping by disabling accounting for parallax correction when generating "
            "resolution masks."
    .type = bool

  resolution_range = None
    .multiple = true
    .type = floats(2)
    .help = "an untrusted resolution range"

  untrusted
    .multiple = True
  {

    panel = None
      .type = int
      .help = "The panel number"
      .help = "If no geometric attributes are given (circle, rectangle, etc)"
      .help = "then the whole panel is masked out"

    circle = None
      .type = ints(3)
      .help = "An untrusted circle (xc, yc, r)"

    rectangle = None
      .type = ints(4)
      .help = "An untrusted rectangle (x0, x1, y0, y1)"

    polygon = None
      .type = ints(value_min=0)
      .help = "The pixel coordinates (fast, slow) that define the corners "
              "of the untrusted polygon. Spots whose centroids fall within "
              "the bounds of the untrusted polygon will be rejected."

    pixel = None
      .type = ints(2, value_min=0)
      .help = "An untrusted pixel (y, x)"

  }

  ice_rings {
    filter = False
      .type = bool
    unit_cell = 4.498,4.498,7.338,90,90,120
      .type = unit_cell
      .help = "The unit cell to generate d_spacings for powder rings."
      .expert_level = 1
    space_group = 194
      .type = space_group
      .help = "The space group used to generate d_spacings for powder rings."
      .expert_level = 1
    width = 0.002
      .type = float(value_min=0.0)
      .help = "The width of an ice ring (in 1/d^2)."
      .expert_level = 1
    d_min = None
      .type = float(value_min=0.0)
      .help = "The high resolution limit (otherwise use detector d_min)"
      .expert_level = 1
  }
""",
    process_includes=True,
)

CacheInfo = namedtuple("CacheInfo", ["hits", "misses", "maxsize", "currsize"])


def lru_equality_cache(maxsize=10):
    """LRU cache that compares keys based on equality... inefficiently.

    Used for dxtbx models that don't have unique id values so can't be
    cached with a normal lru_cache.

    Args:
        maxsize (int): The maximum number of old results to remember
    """

    def _decorator(f):
        class Scope:
            pass

        cache_data = Scope()
        cache_data.cache = []
        cache_data.hits = 0
        cache_data.misses = 0

        def _wrapper_function(*args, **kwargs):
            for i, (key_args, key_kwargs, key_result) in enumerate(cache_data.cache):
                if key_args == args and key_kwargs == kwargs:
                    cache_data.cache.append(cache_data.cache.pop(i))
                    cache_data.hits += 1
                    return key_result
            result = f(*args, **kwargs)
            cache_data.misses += 1
            cache_data.cache.append((args, kwargs, result))
            if len(cache_data.cache) > maxsize:
                cache_data.cache = cache_data.cache[1:]
            return result

        def _generate_cache_info():
            return CacheInfo(
                hits=cache_data.hits,
                misses=cache_data.misses,
                maxsize=maxsize,
                currsize=len(cache_data.cache),
            )

        _wrapper_function.__wrapped__ = f
        _wrapper_function.cache_info = _generate_cache_info

        return _wrapper_function

    return _decorator


def generate_ice_ring_resolution_ranges(beam, panel, params):
    """
    Generate a set of resolution ranges from the ice ring parameters
    """
    if params.filter is True:

        # Get the crystal symmetry
        crystal_symmetry = crystal.symmetry(
            unit_cell=params.unit_cell, space_group=params.space_group.group()
        )

        # Get the half width
        half_width = params.width * 0.5

        # Set the high resolution
        if params.d_min is None:
            d_min = panel.get_max_resolution_at_corners(beam.get_s0())
        else:
            d_min = params.d_min

        # Build the miller set
        ms = crystal_symmetry.build_miller_set(anomalous_flag=False, d_min=d_min)
        ms = ms.sort(by_value="resolution")

        # Yield all the d ranges
        for j, d in enumerate(ms.d_spacings().data()):
            d_sq_inv = 1.0 / (d**2)
            d_sq_inv_min = d_sq_inv - half_width
            d_sq_inv_max = d_sq_inv + half_width
            d_max = math.sqrt(1.0 / d_sq_inv_min)
            d_min = math.sqrt(1.0 / d_sq_inv_max)
            yield (d_min, d_max)


@lru_equality_cache(maxsize=3)
def _get_resolution_masker(beam, panel):
    t0 = time.perf_counter()
    masker = ResolutionMaskGenerator(beam, panel)
    t1 = time.perf_counter()
    logger.debug(f"ResolutionMaskGenerator calculation took {t1 - t0:.4f} seconds")
    return masker


def _apply_resolution_mask(mask, beam, panel, *args):
    _get_resolution_masker(beam, panel).apply(mask, *args)


def generate_mask(
    imageset: ImageSet, params: libtbx.phil.scope_extract
) -> Tuple[flex.bool]:
    """Generate a mask based on the input parameters.

    Args:
      imageset (ImageSet): The imageset for which to generate a mask
      params (libtbx.phil.scope_extract): The phil parameters for mask generation. This
        should be an extract of `dials.util.masking.phil_scope`
    """
    # Get the detector and beam
    detector = imageset.get_detector()
    beam = imageset.get_beam()

    # Create the mask for each panel
    masks = []
    for index, panel in enumerate(detector):

        mask = flex.bool(flex.grid(reversed(panel.get_image_size())), True)

        # Add a border around the image
        if params.border > 0:
            logger.debug(f"Generating border mask:\n border = {params.border}")
            border = params.border
            height, width = mask.all()
            borderx = flex.bool(flex.grid(border, width), False)
            bordery = flex.bool(flex.grid(height, border), False)
            mask[0:border, :] = borderx
            mask[-border:, :] = borderx
            mask[:, 0:border] = bordery
            mask[:, -border:] = bordery

        # Apply the untrusted regions
        for region in params.untrusted:
            if region.panel is None:
                region.panel = 0
            if region.panel == index:
                if not any(
                    [region.circle, region.rectangle, region.polygon, region.pixel]
                ):
                    mask[:, :] = flex.bool(flex.grid(mask.focus()), False)
                    continue

                if region.circle is not None:
                    xc, yc, radius = region.circle
                    logger.debug(
                        "Generating circle mask:\n"
                        + f" panel = {region.panel}\n"
                        + f" xc = {xc}\n"
                        + f" yc = {yc}\n"
                        + f" radius = {radius}"
                    )
                    mask_untrusted_circle(mask, xc, yc, radius)
                if region.rectangle is not None:
                    x0, x1, y0, y1 = region.rectangle
                    logger.debug(
                        "Generating rectangle mask:\n"
                        + f" panel = {region.panel}\n"
                        + f" x0 = {x0}\n"
                        + f" y0 = {y0}\n"
                        + f" x1 = {x1}\n"
                        + f" y1 = {y1}"
                    )
                    mask_untrusted_rectangle(mask, x0, x1, y0, y1)
                if region.polygon is not None:
                    assert (
                        len(region.polygon) % 2 == 0
                    ), "Polygon must contain 2D coords"
                    vertices = []
                    for i in range(int(len(region.polygon) / 2)):
                        x = region.polygon[2 * i]
                        y = region.polygon[2 * i + 1]
                        vertices.append((x, y))
                    polygon = flex.vec2_double(vertices)
                    logger.debug(
                        f"Generating polygon mask:\n panel = {region.panel}\n"
                        + "\n".join(f" coord = {vertex}" for vertex in vertices)
                    )
                    mask_untrusted_polygon(mask, polygon)
                if region.pixel is not None:
                    mask[region.pixel] = False

        # PxMmStrategy to use for generating resolution masks
        if params.disable_parallax_correction:
            panel = copy.deepcopy(panel)
            panel.set_px_mm_strategy(SimplePxMmStrategy())

        # Generate high and low resolution masks
        if params.d_min is not None:
            logger.debug(f"Generating high resolution mask:\n d_min = {params.d_min}")
            _apply_resolution_mask(mask, beam, panel, 0, params.d_min)
        if params.d_max is not None:
            logger.debug(f"Generating low resolution mask:\n d_max = {params.d_max}")
            d_max = params.d_max
            d_inf = max(d_max + 1, 1e9)
            _apply_resolution_mask(mask, beam, panel, d_max, d_inf)

        try:
            # Mask out the resolution range
            for drange in params.resolution_range:
                d_min = min(drange)
                d_max = max(drange)
                assert d_min < d_max, "d_min must be < d_max"
                logger.debug(
                    "Generating resolution range mask:\n"
                    + f" d_min = {d_min}\n"
                    + f" d_max = {d_max}"
                )
                _apply_resolution_mask(mask, beam, panel, d_min, d_max)
        except TypeError:
            # Catch the default value None of params.resolution_range
            if any(params.resolution_range):
                raise

        # Mask out the resolution ranges for the ice rings
        for drange in generate_ice_ring_resolution_ranges(
            beam, panel, params.ice_rings
        ):
            d_min = min(drange)
            d_max = max(drange)
            assert d_min < d_max, "d_min must be < d_max"
            logger.debug(
                "Generating ice ring mask:\n"
                + f" d_min = {d_min:.4f}\n"
                + f" d_max = {d_max:.4f}"
            )
            _apply_resolution_mask(mask, beam, panel, d_min, d_max)

        # Add to the list
        masks.append(mask)

    # Return the mask
    return tuple(masks)

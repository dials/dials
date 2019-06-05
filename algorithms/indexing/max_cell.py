from __future__ import absolute_import, division, print_function

import logging

from scitbx.array_family import flex
from dials.algorithms.spot_finding.per_image_analysis import ice_rings_selection
from dials.algorithms.indexing.nearest_neighbor import NeighborAnalysis
from dials.algorithms.indexing import DialsIndexError

logger = logging.getLogger(__name__)


def find_max_cell(
    reflections,
    max_cell_multiplier=1.3,
    step_size=45,
    nearest_neighbor_percentile=None,
    histogram_binning="linear",
    nn_per_bin=5,
    max_height_fraction=0.25,
    filter_ice=True,
    filter_overlaps=True,
    overlaps_border=0,
):
    logger.debug("Finding suitable max_cell based on %i reflections" % len(reflections))
    # Exclude potential ice-ring spots from nearest neighbour analysis if needed
    if filter_ice:

        ice_sel = ice_rings_selection(reflections)
        reflections = reflections.select(~ice_sel)
        logger.debug(
            "Rejecting %i reflections at ice ring resolution" % ice_sel.count(True)
        )

    # need bounding box in reflections to find overlaps; this is not there if
    # spots are from XDS (for example)
    if filter_overlaps and "bbox" in reflections:
        overlap_sel = flex.bool(len(reflections), False)

        overlaps = reflections.find_overlaps(border=overlaps_border)
        for item in overlaps.edges():
            i0 = overlaps.source(item)
            i1 = overlaps.target(item)
            overlap_sel[i0] = True
            overlap_sel[i1] = True
        logger.debug(
            "Rejecting %i overlapping bounding boxes" % overlap_sel.count(True)
        )
        reflections = reflections.select(~overlap_sel)
    logger.debug("%i reflections remain for max_cell identification" % len(reflections))

    if not len(reflections):
        raise DialsIndexError(
            "Too few spots remaining for nearest neighbour analysis (%d)"
            % len(reflections)
        )

    NN = NeighborAnalysis(
        reflections,
        step_size=step_size,
        max_height_fraction=max_height_fraction,
        tolerance=max_cell_multiplier,
        percentile=nearest_neighbor_percentile,
        histogram_binning=histogram_binning,
        nn_per_bin=nn_per_bin,
    )

    return NN

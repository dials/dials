from __future__ import annotations

import boost_adaptbx.boost.python

ext = boost_adaptbx.boost.python.import_ext("dials_model_data_ext")

# The flex import is required to use the c-grid conversions for the Shoebox mask
# array, which is of type uint8_t. This is currently exposed in
# dials_array_family_flex_ext, though this may be better moved elsewhere
# - potentially to scitbx alongside the other c-grid conversions like int.
# FIXME - to be discussed.
from dials.array_family import flex  # noqa: F401, E402
from dials_model_data_ext import *  # noqa: F403, E402

__all__ = (  # noqa: F405
    "AdjacencyList",
    "AdjacentVerticesIter",
    "Centroid",
    "CentroidData",
    "EdgeDescriptor",
    "ImageDouble",
    "ImageInt",
    "ImageVolume",
    "Intensity",
    "IntensityData",
    "MultiPanelImageVolume",
    "Observation",
    "PixelList",
    "PixelListLabeller",
    "PositionData",
    "Prediction",
    "Ray",
    "Shoebox",
    "make_image",
)

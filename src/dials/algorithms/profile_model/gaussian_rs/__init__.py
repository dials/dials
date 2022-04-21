from __future__ import annotations

import dials.algorithms.profile_model.modeller  # noqa: F401; lgtm; true import dependency
from dials.algorithms.profile_model.gaussian_rs.model import Model, phil_scope
from dials_algorithms_profile_model_gaussian_rs_ext import (
    BBoxCalculator2D,
    BBoxCalculator3D,
    BBoxCalculatorIface,
    BBoxMultiCalculator,
    CoordinateSystem,
    CoordinateSystem2d,
    GaussianRSProfileModeller,
    MaskCalculator2D,
    MaskCalculator3D,
    MaskCalculatorIface,
    MaskMultiCalculator,
    PartialityCalculator2D,
    PartialityCalculator3D,
    PartialityCalculatorIface,
    PartialityMultiCalculator,
    ideal_profile_double,
    ideal_profile_float,
    zeta_factor,
)

__all__ = [
    "BBoxCalculator",
    "BBoxCalculator2D",
    "BBoxCalculator3D",
    "BBoxCalculatorIface",
    "BBoxMultiCalculator",
    "CoordinateSystem",
    "CoordinateSystem2d",
    "GaussianRSProfileModeller",
    "MaskCalculator",
    "MaskCalculator2D",
    "MaskCalculator3D",
    "MaskCalculatorIface",
    "MaskMultiCalculator",
    "Model",
    "PartialityCalculator",
    "PartialityCalculator2D",
    "PartialityCalculator3D",
    "PartialityCalculatorIface",
    "PartialityMultiCalculator",
    "ideal_profile_double",
    "ideal_profile_float",
    "phil_scope",
    "zeta_factor",
]


def BBoxCalculator(crystal, beam, detector, goniometer, scan, delta_b, delta_m):
    """Return the relevant bbox calculator."""
    if goniometer is None or scan is None or scan.is_still():
        algorithm = BBoxCalculator2D(beam, detector, delta_b, delta_m)
    else:
        algorithm = BBoxCalculator3D(beam, detector, goniometer, scan, delta_b, delta_m)
    return algorithm


def PartialityCalculator(crystal, beam, detector, goniometer, scan, sigma_m):
    """Return the relevant partiality calculator."""
    if goniometer is None or scan is None or scan.is_still():
        algorithm = PartialityCalculator2D(beam, sigma_m)
    else:
        algorithm = PartialityCalculator3D(beam, goniometer, scan, sigma_m)
    return algorithm


def MaskCalculator(crystal, beam, detector, goniometer, scan, delta_b, delta_m):
    """Return the relevant partiality calculator."""
    if goniometer is None or scan is None or scan.is_still():
        algorithm = MaskCalculator2D(beam, detector, delta_b, delta_m)
    else:
        algorithm = MaskCalculator3D(beam, detector, goniometer, scan, delta_b, delta_m)
    return algorithm

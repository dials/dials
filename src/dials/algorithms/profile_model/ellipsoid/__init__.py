from __future__ import annotations

from math import pi, sqrt

from dials.array_family import flex  # noqa: F401;
from dials_algorithms_profile_model_ellipsoid_ext import *  # noqa: F401, F403;


def mosaicity_from_eigen_decomposition(eigen_values):
    return tuple(sqrt(e) * 180.0 / pi if e > 0 else 0.0 for e in eigen_values)

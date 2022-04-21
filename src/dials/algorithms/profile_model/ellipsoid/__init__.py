from __future__ import annotations

from math import pi, sqrt

from dials.array_family import flex  # noqa: F401;
from dials_algorithms_profile_model_ellipsoid_ext import *  # noqa: F401, F403;


def mosaicity_from_eigen_decomposition(eigen_values):
    return (
        sqrt(eigen_values[0]) * 180.0 / pi,
        sqrt(eigen_values[1]) * 180.0 / pi,
        sqrt(eigen_values[2]) * 180.0 / pi,
    )

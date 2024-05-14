from __future__ import annotations

# Import dataclass
from dataclasses import dataclass
from typing import Tuple

import matplotlib
import numpy as np

from dials.algorithms.beam_position.helper_functions import (
    Line2D,
    PlotParams,
    normalize,
    plot_profile,
)

matplotlib.use("Agg")


@dataclass
class InversionMethodParams:
    """
    Parameters for the inversion method

    Parameters
    ----------
    guess_position : Tuple[int, int], optional
        The guess beam position (x0, y0) in pixels. Inversion calculation
        is performed only in a window centered around x0 (or y0).
        If not provided, the method will set (x0, y0) to the image center.
    inversion_window_width : int, optional
        The width of the inversion window in pixels. Default is 50 (pixels).
    bad_pixel_threshold : int, optional
        The threshold for bad pixels. All pixels above this value will be set
        to zero. Default is 20000.
    plot : bool, optional
        Whether to plot the inversion results. Default is False.
    verbose : bool, optional
        Whether to print beam position to the output. Default is True.
    filename : str, optional
        The filename for the plot. Default is 'fig.png'.
    """

    guess_position: Tuple[int, int] = None
    inversion_window_width: int = 50
    bad_pixel_threshold: int = 20000
    plot: bool = False
    verbose: bool = True
    filename: str = "fig.png"


def beam_position_from_inversion(
    image: np.ndarray,
    params: InversionMethodParams,
) -> Tuple[float, float]:

    """
    Compute the beam position using the inversion method

    Parameters
    ----------
    image : 2D numpy.ndarray
        The diffraction image.
    params : InversionMethodParams
        Parameters for the inversion method.

    Returns
    -------
    beam_position : Tuple[float, float]
        The beam position (x, y) in pixels.
    """

    image[image > params.bad_pixel_threshold] = 0

    data_x = find_inversion_max(image, params, axis="x")
    data_y = find_inversion_max(image, params, axis="y")

    bx = data_x["beam_position"]
    by = data_y["beam_position"]

    if params.verbose:
        print(f"Beam position from inversion: ({bx:.2f}, {by:.2f})")

    if params.plot:

        i1, i2 = data_x["bin_position"]
        j1, j2 = data_y["bin_position"]
        span_xy = i1, i2, j1, j2

        p = PlotParams(
            image,
            profiles_x=data_x["profile"],
            profiles_y=data_y["profile"],
            beam_position=(bx, by),
            span_xy=span_xy,
            filename=params.filename,
        )

        plot_profile(p)

    return bx, by


def find_inversion_max(image, params: InversionMethodParams, axis="x") -> dict:
    """
    Compute the beam position using the inversion method

    Parameters
    ----------
    image : 2D numpy.ndarray
        The diffraction image.
    params : InversionMethodParams
        Parameters for the inversion method.
    axis : str, optional
        Either 'x' or 'y'.

    Returns
    -------
    data : dict
        A dictionary containing projected profile and beam position.
    """

    if params.guess_position is not None:
        if len(params.guess_position) != 2:
            raise ValueError("Beam position must be a tuple of two ints")

    if params.guess_position is not None:
        if axis == "x":
            center = params.guess_position[0]
        else:
            center = params.guess_position[1]
    else:
        center = None

    if axis == "x":
        profile = image[:, :].max(axis=0)
    elif axis == "y":
        profile = image[:, :].max(axis=1)
    else:
        msg = f"Unknown projection axis '{axis}'. Use either 'x' or 'y'."
        raise ValueError(msg)

    n = len(profile)
    if center is None:
        center = int(n / 2)

    profile = normalize(profile)
    profile_indices = np.arange(n)

    indices = np.arange(
        center - params.inversion_window_widthwidth,
        center + params.inversion_window_width,
        1,
    )
    correlations = []
    for index in indices:
        correlations.append(invert_and_correlate(profile, index))

    correlations = np.array(correlations)
    index = correlations.argmax()

    data = {}

    line_profile = Line2D(profile_indices, profile, c="C0", lw=1.0)
    inversion_profile = Line2D(indices, correlations, c="C3", lw=0.5)
    lines = [line_profile, inversion_profile]

    data["profile"] = lines
    data["beam_position"] = indices[0] + index

    return data


def invert_and_correlate(x, index):
    """Given an 1D array x, compute inverted array inv_x
    (inversion around an element with index 'index')
    and return a sum of (x * inv_x).
    """

    inv_x = np.zeros(len(x))
    n = len(x)
    half = int(n / 2)

    # Compute the inverted 1D array
    if index <= half:
        left = x[0:index]
        right = x[index : 2 * index]
        inv_x[0:index] = right[::-1]
        inv_x[index : 2 * index] = left[::-1]
    else:
        right = x[index:]
        width = len(right)
        left = x[index - width : index]
        inv_x[index - width : index] = right[::-1]
        inv_x[index:] = left[::-1]

    return np.sum(x * inv_x)

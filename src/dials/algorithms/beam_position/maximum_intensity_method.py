from __future__ import annotations

from dataclasses import dataclass
from typing import Tuple

import numpy as np
from plot_params import PlotParams
from plot_profile import plot_profile

from dials.algorithms.beam_position.helper_functions import Line2D, normalize, smooth


@dataclass
class MaxMethodParams:
    """
    Parameters for the max method

    Parameters
    ----------
    bad_pixel_threshold : int, optional
        Set all pixels above this value to zero. Default is 20000.
    bin_step : int, optional
        Distance (in pixels) between neighboring bins,
    bin_width : int, optional
        The width of the bin used to find the region of max intensity (pixels).
    convolution_width : int, optional
        The width of the convolution kernel used for smoothing.
        Default is 2 (pixels).
    n_convolutions : int, optional
        Number of consecutive times the projected profile is smoothed.
        Default is 2.
    plot : bool, optional
        Plot the diffraction image with the computed beam position
        Default is False.
    verbose : bool, optional
        Print out the computed beam position. Default is True.
    filename : str, optional
        Filename to save the plot. Default is 'fig.png'.
    """

    bad_pixel_threshold: int = 20000
    bin_step: int = 10
    bin_width: int = 20
    convolution_width: int = 2
    n_convolutions: int = 2
    plot: bool = False
    verbose: bool = True
    filename: str = "fig.png"


def beam_position_from_max(
    image: np.ndarray,
    params: MaxMethodParams,
) -> Tuple[float, float]:
    """
    Compute beam position using the maximum pixel method

    Parameters
    ----------
    image : 2D numpy.ndarray
        The diffraction image.
    params : MaxMethodParams
        Parameters for the max method.

    Returns
    -------
    bx, by: Tuple[float, float]
        The beam position in the x and y directions.
    """

    image[image > params.bad_pixel_threshold] = 0

    data_x = find_max(image, params, axis="x")
    data_y = find_max(image, params, axis="y")

    bx = data_x["beam_position"]
    by = data_y["beam_position"]

    if params.plot:
        i1, i2 = data_x["bin_position"]
        j1, j2 = data_y["bin_position"]
        span_xy = i1, i2, j1, j2
        p = PlotParams(
            image,
            profiles_x=data_x["profile"],
            profiles_y=data_y["profile"],
            span_xy=span_xy,
            beam_position=(bx, by),
            filename=params.filename,
        )
        plot_profile(p)

    return bx, by


def max_intensity_binning(
    profile_smooth: np.ndarray, profile_max: np.ndarray, params: MaxMethodParams
) -> Tuple[int, int, int]:
    """
    Determine the maximum intensity region using binning

    Parameters
    ----------
    profile_max : 1D numpy.ndarray
        The projected profile of maximum pixels.
    profile_smooth : 1D numpy.ndarray
        The projected average profile after smoothing.
    params : MaxMethodParams
        Parameters for the max method.

    Returns
    -------
    beam_position, i1, i2 : Tuple[float, int, int]
        The beam position and the indices of the maximum intensity region.
    """

    n = len(profile_max)
    n_end = (n // params.bin_width) * params.bin_width
    bins = np.arange(0, n_end, params.bin_step)

    bin_values = []
    bin_indices = []
    for ibin in bins:
        bvalue = profile_smooth[ibin : ibin + params.bin_width]
        bin_indices.append((ibin, ibin + params.bin_width))
        bin_values.append(bvalue.sum())

    bin_values = np.array(bin_values)
    max_index = np.argmax(bin_values)

    i1, i2 = bin_indices[max_index]
    selected = np.array(profile_max)
    selected[0:i1] = 0
    selected[i2:] = 0

    beam_position = np.argmax(selected)

    return beam_position, i1, i2


def find_max(image: np.ndarray, params: MaxMethodParams, axis: str = "x") -> dict:
    """ "
    Project the diffraction image and determine the beam position using the
    maximum pixel method

    Note
    ----
    The method projects a 2D diffraction image into x or y direction
    and convolutes (smooths) the 1D projected profile. Next, it does some
    binning of the 1D data to determine the region with largest intensity.
    The beam position is then determined as the maximal pixel within the
    region of the largest intensity.

    Parameters
    ----------
    image : 2D numpy.ndarray
        The diffraction image.
    params : MaxMethodParams
        Parameters for the max method.
    axis : str, optional
        Either 'x' or 'y' to determine the direction of projection.
        Default is 'x'.

    Returns
    -------
    data : dict
        Data about projected profiles and beam position.
    """

    if axis == "x":
        axis_index = 0
    elif axis == "y":
        axis_index = 1
    else:
        msg = f"Unknown projection axis '{axis}'. Use 'x' or 'y'."
        raise ValueError(msg)

    profile_max = image[:, :].max(axis=axis_index)
    profile_smooth = image[:, :].mean(axis=axis_index)

    for _ in range(params.n_convolutions):
        profile_smooth = smooth(profile_smooth, params.convolution_width)

    profile_max = normalize(profile_max)
    profile_smooth = normalize(profile_smooth)

    line_max = Line2D(x=np.arange(len(profile_max)), y=profile_max, c="C0", lw=1)
    line_smooth = Line2D(
        x=np.arange(len(profile_smooth)), y=profile_smooth, c="C3", lw=0.5
    )

    beam_position, i1, i2 = max_intensity_binning(profile_smooth, profile_max, params)
    data = {}
    data["profiles"] = [line_max, line_smooth]
    data["bin_position"] = i1, i2
    data["beam_position"] = beam_position

    return data

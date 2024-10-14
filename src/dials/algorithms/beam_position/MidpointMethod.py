"""Define a class that searches for beam position using midpoint method"""
from __future__ import annotations

import random
from dataclasses import dataclass
from typing import List, Optional, Tuple

import numpy as np

from dials.algorithms.beam_position.helper_functions import (
    Line2D,
    PlotParams,
    normalize,
    plot_profile,
    remove_percentiles,
    smooth,
)


class MidpointMethodSolver:
    
    def __init__(self, midpoint_range=(0.3, 0.9, 0.01), exclude_range=None):
        """
        Finds beam position from a projected profile using midpoint method

        Parameters
        ----------
        midpoint_range : Tuple[int, int, int], optional
            Midpoint slicing range (start, stop, step). The projected profile
            is normalized to maximal height of one. `start` and `stop` should
            be be between 0 and 1 since data is normalized. `step` should be
            positive.  Default is (start=0.3, stop=0.9, step=0.01).
        exclude_range : List[Tuple[int, int]], optional
            List of pixel ranges of the form (start, stop) to exclude from
            the projected profile. Default is None.
        """

        self.midpoint_range = midpoint_range
        self.exclude_range = exclude_range

    def find_beam_position(self, projected_profile, plot_axes=None):
        """
        Beam position from a projected profile using the midpoint method

        Parameters
        ----------
        projected_profile : 2D numpy.ndarray
            A single diffraction image.
        Returns
        -------
        x, y : Tuple[float, float]
            The beam center position in pixels (x, y).
    """

        # img_clean = remove_percentiles(image, percentile=1 - discard_percentile * 0.01)
        # profile, midpoints, levels_x = find_midpoint(img_clean, params, axis="x")

        y_intersections, groups_of_midpoints = find_midpoints(profile) 
        average_midpoint = pick_by_occurrence(groups_of_midpoints)

        if plot_axes:

            indices_x = np.arange(len(profile_x))
            line_x = [Line2D(indices_x, profile_x)]

#        for midpoints, levels in zip(midpoints_x, levels_x):
#            rand_color = (random.random(), random.random(), random.random())
#            line_x.append(
#                Line2D(midpoints, levels, c=rand_color, lw=0.0, marker="o", ms=0.5)
#            )
#
#        indices_y = np.arange(len(profile_y))
#        line_y = [Line2D(profile_y, indices_y)]
#
#        for midpoints, levels in zip(midpoints_y, levels_y):
#
#            rand_color = (random.random(), random.random(), random.random())
#
#            line_y.append(
#                Line2D(levels, midpoints, c=rand_color, lw=0.0, marker="o", ms=0.5)
#            )
#
#        p = PlotParams(
#            image,
#            profiles_x=line_x,
#            profiles_y=line_y,
#            beam_position=(x0, y0),
#            span_xy=None,
#            filename=plot_filename,
#        )
#
#        plot_profile(p)
#
        return x0, y0



def pick_by_occurrence(midpoint_groups: List[List[float]]) -> float:
    """ 
    Compute the average midpoint value

    Given a list of midpoint groups (a group of midpoints 
    is just a list of midpoints with similar values) find the group with
    maximum number of midpoints and return the average midpoint value for
    that group.
    
    Arguments
    ---------
    midpoint_groups : List[List[float]]
        A list of midpoint groups (a midpoint group is a lists of floats).
    
    Returns
    -------
    average : float 
        The arithmetic average of midpoints from the group with highest
        number of elements.
    
    NOTE: The midpoint groups are already sorted according to some other 
    criterium, so we are only looking at the first three (or less) groups. 
    """

    occurences = [len(group) for group in midpoint_groups]

    if len(midpoint_groups) >= 3:
        max_occurence = max(occurences[0:3])
        max_index = occurences.index(max_occurence)
        average = np.mean(midpoint_groups[max_index])
    elif len(midpoint_groups) == 2:
        max_occurence = max(occurences[0:2])
        max_index = occurences.index(max_occurence)
        average = np.mean(midpoint_groups[max_index])
    else:
        average = np.mean(midpoint_groups[0])
    return average



def add_peak_and_width(
    peaks: List[List[float]],
    widths: List[List[float]],
    levels: List[List[float]],
    m: Tuple[float, float, float],
    threshold=40,
):

    new_peak, new_width, ycut = m

    # Iterate through existing peaks and widths
    for i, peak_list in enumerate(peaks):

        avg_peak = sum(peak_list) / len(peak_list)

        if abs(new_peak - avg_peak) <= threshold:
            # Add new peak to the sublist and corresponding width to widths list
            peak_list.append(new_peak)
            widths[i].append(new_width)
            levels[i].append(ycut)
            return  # Exit function after adding peak and width

    # If no sublist satisfies the condition, create a new sublist
    peaks.append([new_peak])
    widths.append([new_width])
    levels.append([ycut])

    return


def sort_peak_by_occurence(
    peaks: List[List[float]], widths: List[List[float]], levels: List[List[float]]
) -> Tuple[List[List[float]], List[float], List[List[float]]]:
    """
    Sort peaks by average width
    """

    average_widths = [sum(width_list) / len(width_list) for width_list in widths]
    zipped_data = list(zip(peaks, levels, average_widths))

    sorted_data = sorted(zipped_data, key=lambda x: x[2], reverse=True)

    sorted_peaks, sorted_levels, sorted_average_widths = zip(*sorted_data)

    return sorted_peaks, sorted_levels, sorted_average_widths


def find_midpoint(
    image: np.ndarray, params: MidpointMethodParams, axis: str = "x"
) -> Tuple[np.ndarray, float]:
    """
    Project the diffraction image and determine the beam center using the
    midpoint method

    Parameters
    ---------
    image : 2D numpy.ndarray
        A diffraction image.
    params : MidpointMethodParams
        Parameters for the midpoint method.
    axis : str, optional
        Either 'x' or 'y' to get the projected profile. Default is 'x'.

    Returns
    -------
    profile, avg_midpoint: Tuple[np.ndarray, float]
        The `profile` is the averaged and convoluted projection of the
        image data along the specified axis, while the `avg_midpoint`
        is the average midpoint computed in the `data_slice`
        selected range.
    """

    if axis == "x":
        profile = image[:, :].mean(axis=0)
        exclude_range = params.exclude_range_x
    elif axis == "y":
        profile = image[:, :].mean(axis=1)
        exclude_range = params.exclude_range_y
    else:
        msg = f"Unknown projection axis '{axis}'. Use either 'x' or 'y'."
        raise ValueError(msg)

    profile[profile < 0] = 0  # Kill negative pixels

    profile = smooth(profile, width=params.convolution_width)
    profile = normalize(profile)

    start, stop, step = params.data_slice
    levels = np.arange(start, stop, step)

    midpoints = []
    widths = []
    levels_out = []

    for level in levels:

        mid_list = middle(profile, level, exclude_range, params.convolution_width)

        for m in mid_list:

            add_peak_and_width(midpoints, widths, levels_out, m)

    midpoints, levels_out, widths = sort_peak_by_occurence(
        midpoints, widths, levels_out
    )
    return profile, midpoints, levels_out


def middle(a, ycut, exclude_range, smooth_width):
    """Compute all the crossings between a and ycut
       and return the middle of the range of the crossings

    Parameters
    ----------
    a : 1D numpy.ndarray
        The 1D array to search for crossings.
    ycut : float
        The y value at which to search for crossings.
    exclude_range : List[Tuple[int, int]]
        A list of tuples defining the ranges to exclude from the search.
    smooth_width : int
        The width of the smoothing window.

    Returns
    -------
    crossings : List[Tuple[int, int, int]]
        A list containing the middle points of all crossings,
        their widths and the ycut.
    """

    # Mark the crossings
    a[a < 0.001] = 0.001

    b = np.array(a)
    b[b > ycut] = -1

    # Mark the excluded regions
    if exclude_range is not None:
        exclude_range = list(exclude_range)
        n_range = int(len(exclude_range) / 2)
        for i in range(n_range):
            start = int(exclude_range[i])
            end = int(exclude_range[i + 1])
            b[start - smooth_width : end + smooth_width] = -2

    transitions = np.where(np.diff(np.sign(b)))[0] + 1

    crossings = []
    for i in range(0, len(transitions), 2):
        start = transitions[i]
        if i + 1 < len(transitions):
            end = transitions[i + 1]
            good_crossing = not (b[start] == -2 or b[end - 1] == -2)

            if good_crossing:
                midpoint = (start + end) / 2
                width = end - start
                if width > 10:
                    crossings.append((midpoint, width, ycut))

    return crossings

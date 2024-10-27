"""Define a class that searches for beam position using midpoint method"""
from __future__ import annotations

from collections import namedtuple
import numpy as np

from dials.algorithms.beam_position.helper_functions import (
    remove_pixels_by_intensity
)
from dials.algorithms.beam_position.project_profile import project


Midpoint = namedtuple('Midpoint', ['x', 'y', 'width'])


class MidpointMethodSolver:

    def __init__(self, image, params, axis='x'):

        percent = params.projection.midpoint.exclude_intensity_percent
        cw = params.projection.midpoint.convolution_width
        self.axis = axis

        if axis == 'x':
            exclude_range = params.projection.exclude_pixel_range_y
            dead_range = params.projection.midpoint.dead_range_x
            if dead_range == 'None':
                dead_range = None

        elif axis == 'y':
            exclude_range = params.projection.exclude_pixel_range_x
            dead_range = params.projection.midpoint.dead_range_y
            if dead_range == 'None':
                dead_range = None

        else:
            raise ValueError(f"Unknown axis: {axis}")

        clean_image = remove_pixels_by_intensity(image, percent=percent)

        profile, max_value = project(clean_image, axis=axis, method='average',
                                     exclude_range=exclude_range,
                                     convolution_width=cw)

        self.params = params
        self.profile = profile
        self.max_value = max_value
        self.dead_range = dead_range

    def find_beam_position(self):

        m = self.params.projection.midpoint
        intersection_range = m.intersection_range
        convolution_width = m.convolution_width
        ignore_width = m.intersection_min_width

        dead_range = self.dead_range
        profile = np.array(self.profile)

        start, stop, step = check_intersection_param(intersection_range)

        intersection_positions = np.arange(start, stop, step)

        midpoint_groups = []

        for level in intersection_positions:

            midpoints = middle(profile, level,
                               dead_range=dead_range,
                               smooth_width=convolution_width,
                               ignore_width=ignore_width)

            for midpoint in midpoints:
                add_midpoint_to_group(midpoint_groups, midpoint)

        sorted_groups = sort_by_average_width(midpoint_groups)

        beam_position = pick_by_occurrence(sorted_groups)

        self.groups_of_midpoints = sorted_groups
        self.beam_position = beam_position

        return beam_position

    def plot(self, figure):

        indices = np.arange(0, len(self.profile))

        if self.axis == 'x':
            ax = figure.axis_x
            ax.axvline(self.beam_position, c='C3', lw=1)
            ax.plot(indices, self.profile, lw=1, c='gray')
            for midpoint_group in self.groups_of_midpoints:
                x_vals = [m.x for m in midpoint_group]
                y_vals = [m.y for m in midpoint_group]
                ax.plot(x_vals, y_vals, marker='o', ms=1, lw=0)
            ax.text(0.01, 0.95, 'method: midpoint', va='top', ha='left',
                    transform=ax.transAxes, fontsize=8)
            label = 'Imax = %.0f' % self.max_value
            ax.text(0.01, 0.75, label, va='top', ha='left',
                    transform=ax.transAxes, fontsize=8)

        elif self.axis == 'y':
            ax = figure.axis_y
            ax.axhline(self.beam_position, c='C3', lw=1)
            ax.plot(self.profile, indices, lw=1, c='gray')
            for midpoint_group in self.groups_of_midpoints:
                y_vals = [m.x for m in midpoint_group]
                x_vals = [m.y for m in midpoint_group]
                ax.plot(x_vals, y_vals, marker='o', ms=1, lw=0)
            ax.text(0.95, 0.99, 'method: midpoint', va='top', ha='right',
                    transform=ax.transAxes, rotation=-90, fontsize=8)
            label = 'Imax = %.0f' % self.max_value
            ax.text(0.75, 0.99, label, va='top', ha='right',
                    transform=ax.transAxes, rotation=-90, fontsize=8)

        else:
            raise ValueError(f"Unknown axis: {self.axis}")


def pick_by_occurrence(midpoint_groups, nmax=3):
    """
    Check which group of midpoints has the highest number and return
    its average position. Consider only the first nmax groups.
    """

    occurences = [len(group) for group in midpoint_groups[0:nmax]]
    max_occurence = max(occurences)
    max_index = occurences.index(max_occurence)
    selected_group = midpoint_groups[max_index]

    average_position = np.array([m.x for m in selected_group]).mean()

    return average_position


def add_midpoint_to_group(groups_of_midpoints, midpoint,
                          distance_threshold=40):

    for i, group in enumerate(groups_of_midpoints):

        avg_x = np.array([m.x for m in group]).mean()

        if abs(midpoint.x - avg_x) <= distance_threshold:
            group.append(midpoint)
            return

    # If midpoint is too far from existing groups, create a new group
    new_group = [midpoint]
    groups_of_midpoints.append(new_group)

    return


def average_width(midpoint_group):
    total_width = sum([midpoint.width for midpoint in midpoint_group])
    return total_width / len(midpoint_group)


def sort_by_average_width(midpoint_groups):

    widths = [average_width(group) for group in midpoint_groups]

    combined = list(zip(widths, midpoint_groups))
    sorted_combined = sorted(combined, key=lambda x: x[0])
    sorted_widths, sorted_groups_of_midpoints = zip(*sorted_combined)

    # Reverse order
    sorted_widths = sorted_widths[::-1]
    sorted_groups_of_midpoints = sorted_groups_of_midpoints[::-1]

    return sorted_groups_of_midpoints


def middle(profile, level, dead_range, smooth_width, ignore_width):
    """Compute midpoints when level line crosses the profile

    Parameters
    ----------
    profile: 1D numpy array of floats
    level: float
        The y value at which to search for midpoints.
    dead_range: list of ints
        A list defining ranges where the beam is hidden.
        For exampe, [256,289,382,522] would define ranges 256-289 and 382-522.
    smooth_width: int
        The width of the smoothing window.
    ignore_width: int
        Ignore all crossings shorter than this width.

    Returns
    -------
    crossings : List of Midpoint tuples
        A list containing the middle points of all crossings.
    """

    profile[profile < 0.001] = 0.001

    b = np.array(profile)
    b[b > level] = -1

    # Mark the dead regions
    if dead_range:
        dead_range = list(dead_range)
        n = int(len(dead_range) / 2)
        for i in range(n):
            start = int(dead_range[i])
            end = int(dead_range[i + 1])
            b[start - smooth_width:end + smooth_width] = -2

    transitions = np.where(np.diff(np.sign(b)))[0] + 1

    crossings = []
    for i in range(0, len(transitions)-1):
        start = transitions[i]
        end = transitions[i + 1]

        positive_crossing = b[start+1] < 0
        good_crossing = not (b[start] == -2 or b[end - 1] == -2)

        if good_crossing and positive_crossing:
            midpoint_position = (start + end) / 2
            width = end - start
            if width > ignore_width:
                point = Midpoint(midpoint_position, level, width)
                crossings.append(point)

    return crossings


def check_intersection_param(midpoint_range):

    if len(midpoint_range) != 3:
        msg = "Midpoint method error! Intersection range requires three "
        msg += " floats (start, stop, step) ranging from 0 to 1"
        raise ValueError(msg)

    start, stop, step = midpoint_range

    for val in midpoint_range:
        if not (isinstance(val, float) or isinstance(val, int)):
            msg = f"The value {val} in intersection range in midpoint method "
            msg + "is neither a float nor an int"
            raise ValueError(msg)

    if (start < 0) or (start > 1):
        msg = "Midpoint method error!\n"
        msg += "The start of the intersection range outside of "
        msg += " the (0, 1) interval."
        raise ValueError(msg)

    if (stop < 0) or (stop > 1):
        msg = "Midpoint method error!\n"
        msg += "The end of the intersection range outside of "
        msg += " the (0, 1) interval."
        raise ValueError(msg)

    return start, stop, step

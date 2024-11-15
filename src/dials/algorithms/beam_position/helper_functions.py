"""Helper functions used to determine the beam position of the imported data"""

from __future__ import annotations

import matplotlib
import numpy as np

matplotlib.use('Agg')


def remove_pixels_by_intensity(image, percent=0.0):

    if percent < 0 or percent >= 100:
        raise ValueError('Exclude intensity percent outside of 0 to 100 range')

    keep_percent = 1.0 - percent * 0.01

    pixels_1d = np.sort(image.flatten())
    ntot = len(pixels_1d)
    ncut = int(ntot * keep_percent)
    image_copy = np.array(image)

    if ncut == ntot:
        return image_copy

    icut = pixels_1d[ncut]
    image_copy[image > icut] = 0
    return image_copy


def normalize(array):
    """Apply a pedestal and normalize a 1D numpy array"""

    min_value = array.min()

    if min_value < 0:
        positive_array = array + abs(min_value)   # Add pedestal
    else:
        positive_array = np.array(array)

    max_value = positive_array.max()

    return positive_array / max_value


def smooth(curve, width=2):
    """
    Smooth a 1D numpy array `curve` with a rectangle convolution.
    The rectangle width is 2*half_width.
    """

    smooth_curve = 0 * curve
    n = len(curve)

    half_width = int(width / 2)

    if half_width <= 0:
        return np.array(curve)

    for i in range(n):

        if i < half_width:
            smooth_curve[i] = curve[0:i + half_width].mean()
        elif i > n - half_width:
            smooth_curve[i] = curve[i - half_width:].mean()
        else:
            smooth_curve[i] = curve[i - half_width:i + half_width].mean()

    return smooth_curve


def get_indices_from_slices(nmax, slices):
    """
    Takes a string of comma-separated numpy slices (e.g. "::2, 2, 3:5:7")
    and turns it into a list of indices. The indices are selected from a range
    [0, nmax] using each slice, and each comma-separated component is included
    in the final array.
    """

    indices = np.array([], dtype=np.dtype('uint'))
    full_range = np.arange(0, nmax, dtype=np.dtype('uint'))

    for slice_str in slices.split(','):
        slice_str = slice_str.replace(' ', '')
        parsed_slice = parse_numpy_slice(slice_str)
        selected = full_range[parsed_slice]
        indices = np.append(indices, selected)

    return np.unique(indices)


def parse_numpy_slice(slice_str):

    if slice_str is None:
        return None
    try:
        return eval(f"np.s_[{slice_str}]")
    except Exception as e:
        raise ValueError(f"Invalid slice format: {slice_str}") from e

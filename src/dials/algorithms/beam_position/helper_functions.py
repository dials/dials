"""Helper functions used to determine the beam position of the imported data"""

from __future__ import annotations

import matplotlib
import numpy as np

matplotlib.use("Agg")


def remove_pixels_by_intensity(image, percent=0.0):
    if percent < 0 or percent >= 100:
        raise ValueError("Exclude intensity percent outside of 0 to 100 range")

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
        positive_array = array + abs(min_value)  # Add pedestal
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
            smooth_curve[i] = curve[0 : i + half_width].mean()
        elif i > n - half_width:
            smooth_curve[i] = curve[i - half_width :].mean()
        else:
            smooth_curve[i] = curve[i - half_width : i + half_width].mean()

    return smooth_curve


def parse_numpy_slice(slice_str):
    if slice_str is None:
        return None
    try:
        return eval(f"np.s_[{slice_str}]")
    except Exception as e:
        raise ValueError(f"Invalid slice format: {slice_str}") from e


def print_progress(image_index, n_images, set_index, n_sets, bar_length=40):
    image_index += 1
    set_index += 1

    percent_images = 1.0 * image_index / n_images
    percent_sets = 1.0 * set_index / n_sets

    move_up = "\033[F"

    img_bar_full = int(percent_images * bar_length)
    image_bar = "=" * img_bar_full + " " * (bar_length - img_bar_full)

    set_bar_full = int(percent_sets * bar_length)
    set_bar = "=" * set_bar_full + " " * (bar_length - set_bar_full)

    bar = f" Set:   [{set_bar}] {100 * percent_sets:0.2f} % "
    bar += f"{set_index:4d}/{n_sets:d}\n"
    bar += f" Image: [{image_bar}] {100 * percent_images:0.2f} % "
    bar += f"{image_index:4d}/{n_images}\n"

    condition_01 = abs(percent_sets - 1.0) < 1.0e-15
    condition_02 = abs(percent_images - 1.0) < 1.0e-15

    if condition_01 and condition_02:
        end_str = ""
    else:
        end_str = f"{move_up}{move_up}\r"

    bar += end_str
    print(bar, end="")

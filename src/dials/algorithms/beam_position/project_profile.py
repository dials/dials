from __future__ import annotations

from collections import namedtuple

import numpy as np

from dials.algorithms.beam_position.helper_functions import normalize as norm
from dials.algorithms.beam_position.helper_functions import smooth

Span = namedtuple('Span', ['min', 'max'])


def project(image, axis="x", method="max", convolution_width=1,
            exclude_range=None, n_convolutions=1, normalize=True):
    """
    Compute image projection on an axis

    Parameters
    ----------
    image: numpy 2D array of floats or ints
        A diffraction image.

    """

    if axis == "x":
        proj_axis = 0
        clean_image = exclude_range_from_image(image, exclude_range, axis="y")
    elif axis == "y":
        proj_axis = 1
        clean_image = exclude_range_from_image(image, exclude_range, axis="x")
    else:
        msg = f"Unknown projection axis '{axis}'. Use either 'x' or 'y'"
        raise ValueError(msg)

    # Still need to include exclude range
    if method == "max":
        profile = clean_image[:, :].max(axis=proj_axis)
    elif method == "average":
        profile = clean_image[:, :].mean(axis=proj_axis)

    for i in range(n_convolutions):
        profile = smooth(profile, width=convolution_width)

    max_value = profile.max()
    min_value = profile.min()
    if normalize:
        profile = norm(profile)

    return profile, max_value, min_value


def convert_range_into_spans(exclude_range):

    spans = []

    if len(list(exclude_range)) == 0:
        return spans

    exclude_range = list(exclude_range[0])

    n_ranges = int(len(exclude_range) / 2)

    mins = exclude_range[0:2*n_ranges]
    maxs = exclude_range[1:2*n_ranges]

    for imin, imax in zip(mins, maxs):

        if imin >= imax:
            raise ValueError(
                "Error in exclude range! Range minimum larger than range "
                f"maximum: {imin} > {imax}"
            )

        spans.append(Span(imin, imax))
    return spans


def exclude_range_from_image(image, exclude_range, axis="x"):

    ny, nx = image.shape
    clean_image = np.array(image)
    spans = convert_range_into_spans(exclude_range)

    for span in spans:
        if span.min < 0 or span.max < 0:
            raise ValueError("Error! Exclude range limits must be positive.")

        if axis == "x" and (span.min > nx or span.max > nx):
            raise ValueError(
                "Error! Exclude range limit larger than image size along x."
            )

        if axis == "y" and (span.min > ny or span.max > ny):
            raise ValueError(
                "Error! Exclude range limit larger than image size along y."
            )

        if axis == "x":
            clean_image[:, span.min:span.max] = 0
        elif axis == 'y':
            clean_image[span.min:span.max, :] = 0

    return clean_image

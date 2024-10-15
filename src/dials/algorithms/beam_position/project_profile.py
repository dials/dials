from dials.algorithms.beam_position.helper_functions import (
    normalize,
    smooth
)


def project(image, axis='x', method='max', exclude_range=None,
            convolution_width=1, n_convolutions=1):

    if axis == "x":
        proj_axis = 0
    elif axis == 'y':
        proj_axis = 1
    else:
        msg = f"Unknown projection axis '{axis}'. Use either 'x' or 'y'"
        raise ValueError(msg)

    # Still need to include exclude range
    if method == 'max':
        profile = image[:, :].max(axis=proj_axis)
    elif method == 'average':
        profile = image[:, :].mean(axis=proj_axis)

    for i in range(n_convolutions):
        profile = smooth(profile, width=convolution_width)

    max_value = profile.max()
    profile = normalize(profile)

    return profile, max_value

from __future__ import absolute_import, division, print_function

import pytest


def generate_image(xsize, ysize):
    from scitbx.array_family import flex

    image = flex.random_double(xsize * ysize)
    image.reshape(flex.grid(ysize, xsize))
    return image


def generate_mask(xsize, ysize):
    from scitbx.array_family import flex

    mask = flex.random_bool(xsize * ysize, 0.9)
    mask.reshape(flex.grid(ysize, xsize))
    return mask


def test_filter():
    from numpy import median
    from dials.algorithms.image.filter import median_filter

    xsize = 200
    ysize = 300
    kernel = (3, 3)

    image = generate_image(xsize, ysize)

    result = median_filter(image, kernel)

    eps = 1e-7

    for j in range(kernel[0], ysize - kernel[0]):
        for i in range(kernel[1], xsize - kernel[1]):
            j0 = j - kernel[0]
            j1 = j + kernel[0] + 1
            i0 = i - kernel[1]
            i1 = i + kernel[1] + 1
            pixels = image[j0:j1, i0:i1]
            value = median(pixels.as_numpy_array())
            assert result[j, i] == pytest.approx(value, abs=eps)


def test_masked_filter():
    from numpy import median
    from dials.algorithms.image.filter import median_filter

    xsize = 200
    ysize = 300
    kernel = (3, 3)

    image = generate_image(xsize, ysize)
    mask = generate_mask(xsize, ysize)
    result = median_filter(image, mask, kernel)
    eps = 1e-7

    for j in range(kernel[0], ysize - kernel[0]):
        for i in range(kernel[1], xsize - kernel[1]):
            if mask[j, i]:
                j0 = j - kernel[0]
                j1 = j + kernel[0] + 1
                i0 = i - kernel[1]
                i1 = i + kernel[1] + 1
                pixels = image[j0:j1, i0:i1]
                pmask = mask[j0:j1, i0:i1]
                pixels = pixels.as_1d().select(pmask.as_1d())
                if len(pixels) & 1:
                    value = median(pixels.as_numpy_array())
                else:
                    pixels = sorted(list(pixels))
                    value = pixels[len(pixels) // 2]
                assert result[j, i] == pytest.approx(value, abs=eps)

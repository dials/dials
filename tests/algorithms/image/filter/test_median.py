from __future__ import annotations

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
                    pixels = sorted(pixels)
                    value = pixels[len(pixels) // 2]
                assert result[j, i] == pytest.approx(value, abs=eps)


@pytest.mark.parametrize("periodic", [False, True])
def test_masked_filter_edges(periodic):
    from dials.algorithms.image.filter import median_filter

    # Sizes that are not powers of 2, so that wrapping negative indices with
    # unsigned arithmetic would pick the wrong pixels
    xsize = 23
    ysize = 37
    kernel = (3, 2)

    image = generate_image(xsize, ysize)
    mask = generate_mask(xsize, ysize)
    result = median_filter(image, mask, kernel, periodic=periodic)

    for j in range(ysize):
        for i in range(xsize):
            pixels = []
            for jj in range(j - kernel[0], j + kernel[0] + 1):
                for ii in range(i - kernel[1], i + kernel[1] + 1):
                    if periodic:
                        jj, ii = jj % ysize, ii % xsize
                    elif not (0 <= jj < ysize and 0 <= ii < xsize):
                        continue
                    if mask[jj, ii]:
                        pixels.append(image[jj, ii])
            expected = sorted(pixels)[len(pixels) // 2] if pixels else 0
            assert result[j, i] == expected

from __future__ import annotations

from math import exp
from random import randint

import pytest
from numpy.random import poisson

from scitbx.array_family import flex

from dials.algorithms.image.threshold import (
    DispersionExtendedThreshold,
    DispersionExtendedThresholdDebug,
    DispersionThreshold,
    DispersionThresholdDebug,
)


class Test:
    def setup_class(self):
        from scitbx.array_family import flex

        spot = flex.double(flex.grid(11, 11))
        for j in range(11):
            for i in range(11):
                spot[j, i] = exp(-((j - 5) ** 2 + (i - 5) ** 2) / 2**2)

        self.image = flex.double(flex.grid(2000, 2000))
        for n in range(200):
            x = randint(0, 2000)
            y = randint(0, 2000)
            for j in range(0, 11):
                for i in range(0, 11):
                    xx = x + i - 5
                    yy = y + j - 5
                    if xx >= 0 and yy >= 0 and xx < 2000 and yy < 2000:
                        self.image[yy, xx] = poisson(100 * spot[j, i])

        background = flex.double(list(poisson(5, 2000 * 2000)))
        background.reshape(flex.grid(2000, 2000))

        self.image += background

        # Create an image
        self.mask = flex.random_bool(2000 * 2000, 0.99)
        self.mask.reshape(flex.grid(2000, 2000))
        self.gain = flex.random_double(2000 * 2000) + 1.0
        self.gain.reshape(flex.grid(2000, 2000))
        self.size = (3, 3)
        self.min_count = 2

    def test_niblack(self):
        from dials.algorithms.image.threshold import niblack

        n_sigma = 3
        niblack(self.image, self.size, n_sigma)

    def test_sauvola(self):
        from dials.algorithms.image.threshold import sauvola

        k = 3
        r = 3
        sauvola(self.image, self.size, k, r)

    def test_index_of_dispersion(self):
        from dials.algorithms.image.threshold import index_of_dispersion

        n_sigma = 3
        index_of_dispersion(self.image, self.size, n_sigma)

    def test_index_of_dispersion_masked(self):
        from dials.algorithms.image.threshold import index_of_dispersion_masked

        n_sigma = 3
        index_of_dispersion_masked(
            self.image, self.mask, self.size, self.min_count, n_sigma
        )

    def test_gain(self):
        from dials.algorithms.image.threshold import gain

        n_sigma = 3
        gain(self.image, self.mask, self.gain, self.size, self.min_count, n_sigma)

    def test_dispersion(self):
        from dials.algorithms.image.threshold import dispersion

        nsig_b = 3
        nsig_s = 3
        dispersion(self.image, self.mask, self.size, nsig_b, nsig_s, self.min_count)

    def test_dispersion_w_gain(self):
        from dials.algorithms.image.threshold import dispersion_w_gain

        nsig_b = 3
        nsig_s = 3
        result1 = dispersion_w_gain(
            self.image, self.mask, self.gain, self.size, nsig_b, nsig_s, self.min_count
        )

        # scaling both the image and the gain should not affect the result
        result2 = dispersion_w_gain(
            2.0 * self.image,
            self.mask,
            2.0 * self.gain,
            self.size,
            nsig_b,
            nsig_s,
            self.min_count,
        )
        assert result1 == result2

        # should get the same result as dispersion if the gain is unity
        from dials.algorithms.image.threshold import dispersion

        result3 = dispersion(
            self.image, self.mask, self.size, nsig_b, nsig_s, self.min_count
        )
        result4 = dispersion_w_gain(
            self.image,
            self.mask,
            (0 * self.gain + 1),
            self.size,
            nsig_b,
            nsig_s,
            self.min_count,
        )
        assert result3 == result4

    def test_dispersion_debug(self):
        from dials.algorithms.image.threshold import (
            DispersionThresholdDebug,
            dispersion,
            dispersion_w_gain,
        )

        nsig_b = 3
        nsig_s = 3
        result1 = dispersion(
            self.image, self.mask, self.size, nsig_b, nsig_s, self.min_count
        )
        debug = DispersionThresholdDebug(
            self.image, self.mask, self.size, nsig_b, nsig_s, 0, self.min_count
        )
        result2 = debug.final_mask()
        assert result1.all_eq(result2)

        result3 = dispersion_w_gain(
            self.image, self.mask, self.gain, self.size, nsig_b, nsig_s, self.min_count
        )
        debug = DispersionThresholdDebug(
            self.image,
            self.mask,
            self.gain,
            self.size,
            nsig_b,
            nsig_s,
            0,
            self.min_count,
        )
        result4 = debug.final_mask()
        assert result3 == result4

    def test_dispersion_threshold(self):
        from dials.algorithms.image.threshold import (
            DispersionThreshold,
            dispersion,
            dispersion_w_gain,
        )
        from dials.array_family import flex

        nsig_b = 3
        nsig_s = 3
        algorithm = DispersionThreshold(
            self.image.all(), self.size, nsig_b, nsig_s, 0, self.min_count
        )

        result1 = dispersion(
            self.image, self.mask, self.size, nsig_b, nsig_s, self.min_count
        )

        result2 = dispersion_w_gain(
            self.image, self.mask, self.gain, self.size, nsig_b, nsig_s, self.min_count
        )

        result3 = flex.bool(flex.grid(self.image.all()))
        result4 = flex.bool(flex.grid(self.image.all()))
        algorithm(self.image, self.mask, result3)
        algorithm(self.image, self.mask, self.gain, result4)
        assert result1 == result3
        assert result2 == result4

    def test_dispersion_extended_threshold(self):
        from dials.algorithms.image.threshold import (
            DispersionExtendedThreshold,
            DispersionExtendedThresholdDebug,
        )
        from dials.array_family import flex

        nsig_b = 3
        nsig_s = 3
        algorithm = DispersionExtendedThreshold(
            self.image.all(), self.size, nsig_b, nsig_s, 0, self.min_count
        )
        result1 = flex.bool(flex.grid(self.image.all()))
        result2 = flex.bool(flex.grid(self.image.all()))
        algorithm(self.image, self.mask, result1)
        algorithm(self.image, self.mask, self.gain, result2)

        debug = DispersionExtendedThresholdDebug(
            self.image, self.mask, self.size, nsig_b, nsig_s, 0, self.min_count
        )
        result3 = debug.final_mask()

        assert result1.all_eq(result3)

        debug = DispersionExtendedThresholdDebug(
            self.image,
            self.mask,
            self.gain,
            self.size,
            nsig_b,
            nsig_s,
            0,
            self.min_count,
        )
        result4 = debug.final_mask()
        assert result2 == result4

    @pytest.mark.parametrize(
        "algorithm", [DispersionThreshold, DispersionExtendedThreshold]
    )
    def test_dispersion_algorithm_symmetry(self, algorithm):

        nsig_b = 3
        nsig_s = 3

        thresholder = algorithm(
            self.image.all(), (3, 3), nsig_b, nsig_s, 0, self.min_count
        )
        result1 = flex.bool(flex.grid(self.image.all()))
        thresholder(self.image, self.mask, result1)

        # Ugh. Is there an easier way to transpose a 2d flex.bool?
        def transpose_a_flex_bool(arr):
            arr_t = arr.as_1d().as_int()
            arr_t.reshape(flex.grid(arr.all()))
            arr_t.matrix_transpose_in_place()
            new_dim = arr_t.all()
            arr_t = arr_t.as_1d().as_bool()
            arr_t.reshape(flex.grid(new_dim))
            return arr_t

        image_t = self.image.matrix_transpose()
        mask_t = transpose_a_flex_bool(self.mask)
        thresholder = algorithm(
            image_t.all(), (3, 3), nsig_b, nsig_s, 0, self.min_count
        )
        result2 = flex.bool(flex.grid(image_t.all()))
        thresholder(image_t, mask_t, result2)
        result2_t = transpose_a_flex_bool(result2)

        assert (result1 == result2_t).all_eq(True)

    @pytest.mark.parametrize(
        "algorithm", [DispersionThresholdDebug, DispersionExtendedThresholdDebug]
    )
    def test_dispersion_debug_algorithm_symmetry(self, algorithm):

        nsig_b = 3
        nsig_s = 3

        threshold = algorithm(
            self.image, self.mask, (3, 3), nsig_b, nsig_s, 0, self.min_count
        )
        result1 = threshold.final_mask()

        # Ugh. Is there an easier way to transpose a 2d flex.bool?
        def transpose_a_flex_bool(arr):
            arr_t = arr.as_1d().as_int()
            arr_t.reshape(flex.grid(arr.all()))
            arr_t.matrix_transpose_in_place()
            new_dim = arr_t.all()
            arr_t = arr_t.as_1d().as_bool()
            arr_t.reshape(flex.grid(new_dim))
            return arr_t

        image_t = self.image.matrix_transpose()
        mask_t = transpose_a_flex_bool(self.mask)
        threshold = algorithm(
            image_t, mask_t, (3, 3), nsig_b, nsig_s, 0, self.min_count
        )
        result2 = threshold.final_mask()
        result2_t = transpose_a_flex_bool(result2)

        assert (result1 == result2_t).all_eq(True)

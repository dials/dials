from __future__ import absolute_import, division, print_function


class Test_Centroid(object):
    EPS = 1e-7

    def setup_class(self):
        self.generate_data()
        self.calculate_gold()

    def test_centroid_points2d(self):
        from dials.algorithms.image.centroid import centroid_points
        from scitbx import matrix

        centroid = centroid_points(self.pixels2d.as_1d(), self.points2d.as_1d())

        assert abs(self.gold2d - matrix.col(centroid.mean())) < self.EPS
        assert abs(self.gold2dvar - matrix.col(centroid.variance())) < self.EPS
        assert (
            abs(self.gold2dubvar - matrix.col(centroid.unbiased_variance())) < self.EPS
        )

    def test_centroid_points3d(self):
        from dials.algorithms.image.centroid import centroid_points
        from scitbx import matrix

        centroid = centroid_points(self.pixels3d.as_1d(), self.points3d.as_1d())

        assert abs(self.gold3d - matrix.col(centroid.mean())) < self.EPS
        assert abs(self.gold3dvar - matrix.col(centroid.variance())) < self.EPS
        assert (
            abs(self.gold3dubvar - matrix.col(centroid.unbiased_variance())) < self.EPS
        )

    def test_centroid_image2d(self):
        from dials.algorithms.image.centroid import centroid_image
        from scitbx import matrix

        centroid = centroid_image(self.pixels2d)

        assert abs(self.gold2d - matrix.col(centroid.mean())) < self.EPS
        assert abs(self.gold2dvar - matrix.col(centroid.variance())) < self.EPS
        assert (
            abs(self.gold2dubvar - matrix.col(centroid.unbiased_variance())) < self.EPS
        )

    def test_centroid_image3d(self):
        from dials.algorithms.image.centroid import centroid_image
        from scitbx import matrix

        centroid = centroid_image(self.pixels3d)

        assert abs(self.gold3d - matrix.col(centroid.mean())) < self.EPS
        assert abs(self.gold3dvar - matrix.col(centroid.variance())) < self.EPS
        assert (
            abs(self.gold3dubvar - matrix.col(centroid.unbiased_variance())) < self.EPS
        )

    def test_centroid_masked_image2d(self):
        from dials.algorithms.image.centroid import centroid_image
        from scitbx import matrix

        centroid = centroid_image(self.pixels2d, self.mask2d)

        assert abs(self.goldmasked2d - matrix.col(centroid.mean())) < self.EPS
        assert abs(self.goldmasked2dvar - matrix.col(centroid.variance())) < self.EPS
        assert (
            abs(self.goldmasked2dubvar - matrix.col(centroid.unbiased_variance()))
            < self.EPS
        )

    def test_centroid_masked_image3d(self):
        from dials.algorithms.image.centroid import centroid_image
        from scitbx import matrix

        centroid = centroid_image(self.pixels3d, self.mask3d)

        assert abs(self.goldmasked3d - matrix.col(centroid.mean())) < self.EPS
        assert abs(self.goldmasked3dvar - matrix.col(centroid.variance())) < self.EPS
        assert (
            abs(self.goldmasked3dubvar - matrix.col(centroid.unbiased_variance()))
            < self.EPS
        )

    def test_centroid_bias(self):
        from dials.algorithms.image.centroid import centroid_image
        from scitbx.array_family import flex

        pixels = flex.double(flex.grid(5, 5), 0)
        pixels[1, 2] = 5
        pixels[2, 2] = 10
        pixels[3, 2] = 5
        pixels[2, 1] = 5
        pixels[2, 3] = 5
        centroid = centroid_image(pixels)
        assert centroid.average_bias_estimate()[0] < 1e-7
        assert centroid.average_bias_estimate()[1] < 1e-7

    @classmethod
    def generate_data(cls):
        from scitbx.array_family import flex
        from random import random, randint

        # Generate a 3d array of pixels and points
        cls.points3d = flex.vec3_double(flex.grid(5, 5, 5))
        cls.pixels3d = flex.double(flex.grid(5, 5, 5))
        cls.mask3d = flex.bool(flex.grid(5, 5, 5))
        for k in range(0, 5):
            for j in range(0, 5):
                for i in range(0, 5):
                    cls.points3d[k, j, i] = (i + 0.5, j + 0.5, k + 0.5)
                    cls.pixels3d[k, j, i] = random()
                    cls.mask3d[k, j, i] = bool(randint(0, 1))

        cls.points2d = flex.vec2_double(flex.grid(5, 5))
        cls.pixels2d = flex.double(flex.grid(5, 5))
        cls.mask2d = flex.bool(flex.grid(5, 5))
        for j in range(0, 5):
            for i in range(0, 5):
                cls.points2d[j, i] = cls.points3d[0, j, i][0:2]
                cls.pixels2d[j, i] = cls.pixels3d[0, j, i]
                cls.mask2d[j, i] = cls.mask3d[0, j, i]

    @classmethod
    def calculate_gold(cls):
        cls.calculate_gold2d()
        cls.calculate_gold3d()
        cls.calculate_gold_masked2d()
        cls.calculate_gold_masked3d()

    @classmethod
    def calculate_gold2d(cls):
        from scitbx.array_family import flex
        from scitbx import matrix

        r_tot = 0.0
        c_tot = 0.0
        d_tot = 0.0

        for (r, c), d in zip(cls.points2d, cls.pixels2d):
            r_tot += d * r
            c_tot += d * c
            d_tot += d

        cls.gold2d = matrix.col((r_tot / d_tot, c_tot / d_tot))
        _r, _c = cls.gold2d

        r_tot = 0.0
        c_tot = 0.0

        for (r, c), d in zip(cls.points2d, cls.pixels2d):
            r_tot += d * (r - _r) ** 2
            c_tot += d * (c - _c) ** 2

        _sr = r_tot / d_tot
        _sc = c_tot / d_tot

        cls.gold2dvar = matrix.col((_sr, _sc))

        # r_tot = 0.0
        # c_tot = 0.0

        # for (r, c), d in zip(cls.points2d, cls.pixels2d):
        #  r_tot += d * (r - _r) ** 2
        #  c_tot += d * (c - _c) ** 2

        # _sr = r_tot / (d_tot-1)
        # _sc = c_tot / (d_tot-1)

        # cls.gold2dubvar = matrix.col((_sr, _sc))

        pixel_x, pixel_y = zip(*cls.points2d)
        xc = flex.mean_and_variance(flex.double(pixel_x), cls.pixels2d.as_1d())
        yc = flex.mean_and_variance(flex.double(pixel_y), cls.pixels2d.as_1d())
        cls.gold2dubvar = matrix.col(
            (xc.gsl_stats_wvariance(), yc.gsl_stats_wvariance())
        )

    @classmethod
    def calculate_gold3d(cls):
        from scitbx.array_family import flex
        from scitbx import matrix

        f_tot = 0.0
        r_tot = 0.0
        c_tot = 0.0
        d_tot = 0.0

        for (f, r, c), d in zip(cls.points3d, cls.pixels3d):
            f_tot += d * f
            r_tot += d * r
            c_tot += d * c
            d_tot += d

        cls.gold3d = matrix.col((f_tot / d_tot, r_tot / d_tot, c_tot / d_tot))

        _f, _r, _c = cls.gold3d

        f_tot = 0.0
        r_tot = 0.0
        c_tot = 0.0

        for (f, r, c), d in zip(cls.points3d, cls.pixels3d):
            f_tot += d * (f - _f) ** 2
            r_tot += d * (r - _r) ** 2
            c_tot += d * (c - _c) ** 2

        _sf = f_tot / d_tot
        _sr = r_tot / d_tot
        _sc = c_tot / d_tot

        cls.gold3dvar = matrix.col((_sf, _sr, _sc))

        # f_tot = 0.0
        # r_tot = 0.0
        # c_tot = 0.0

        # for (f, r, c), d in zip(cls.points3d, cls.pixels3d):
        #  f_tot += d * (f - _f) ** 2
        #  r_tot += d * (r - _r) ** 2
        #  c_tot += d * (c - _c) ** 2

        # _sf = f_tot / (d_tot-1)
        # _sr = r_tot / (d_tot-1)
        # _sc = c_tot / (d_tot-1)

        # cls.gold3dubvar = matrix.col((_sf, _sr, _sc))

        pixel_x, pixel_y, pixel_z = zip(*cls.points3d)
        xc = flex.mean_and_variance(flex.double(pixel_x), cls.pixels3d.as_1d())
        yc = flex.mean_and_variance(flex.double(pixel_y), cls.pixels3d.as_1d())
        zc = flex.mean_and_variance(flex.double(pixel_z), cls.pixels3d.as_1d())
        cls.gold3dubvar = matrix.col(
            (
                xc.gsl_stats_wvariance(),
                yc.gsl_stats_wvariance(),
                zc.gsl_stats_wvariance(),
            )
        )

    @classmethod
    def calculate_gold_masked2d(cls):
        from scitbx.array_family import flex
        from scitbx import matrix

        r_tot = 0.0
        c_tot = 0.0
        d_tot = 0.0

        for (r, c), d, m in zip(cls.points2d, cls.pixels2d, cls.mask2d):
            if m:
                r_tot += d * r
                c_tot += d * c
                d_tot += d

        cls.goldmasked2d = matrix.col((r_tot / d_tot, c_tot / d_tot))

        _r, _c = cls.goldmasked2d

        r_tot = 0.0
        c_tot = 0.0

        for (r, c), d, m in zip(cls.points2d, cls.pixels2d, cls.mask2d):
            if m:
                r_tot += d * (r - _r) ** 2
                c_tot += d * (c - _c) ** 2

        _sr = r_tot / d_tot
        _sc = c_tot / d_tot

        cls.goldmasked2dvar = matrix.col((_sr, _sc))

        # r_tot = 0.0
        # c_tot = 0.0

        # for (r, c), d, m in zip(cls.points2d, cls.pixels2d, cls.mask2d):
        #  if m:
        #    r_tot += d * (r - _r) ** 2
        #    c_tot += d * (c - _c) ** 2

        # _sr = r_tot / (d_tot-1)
        # _sc = c_tot / (d_tot-1)

        # cls.goldmasked2dubvar = matrix.col((_sr, _sc))

        pixel_x = flex.double()
        pixel_y = flex.double()
        pixel_d = flex.double()
        for (x, y), d, m in zip(cls.points2d, cls.pixels2d, cls.mask2d):
            if m:
                pixel_x.append(x)
                pixel_y.append(y)
                pixel_d.append(d)

        xc = flex.mean_and_variance(flex.double(pixel_x), pixel_d)
        yc = flex.mean_and_variance(flex.double(pixel_y), pixel_d)
        cls.goldmasked2dubvar = matrix.col(
            (xc.gsl_stats_wvariance(), yc.gsl_stats_wvariance())
        )

    @classmethod
    def calculate_gold_masked3d(cls):
        from scitbx.array_family import flex
        from scitbx import matrix

        f_tot = 0.0
        r_tot = 0.0
        c_tot = 0.0
        d_tot = 0.0

        for (f, r, c), d, m in zip(cls.points3d, cls.pixels3d, cls.mask3d):
            if m:
                f_tot += d * f
                r_tot += d * r
                c_tot += d * c
                d_tot += d

        cls.goldmasked3d = matrix.col((f_tot / d_tot, r_tot / d_tot, c_tot / d_tot))

        _f, _r, _c = cls.goldmasked3d

        f_tot = 0.0
        r_tot = 0.0
        c_tot = 0.0

        for (f, r, c), d, m in zip(cls.points3d, cls.pixels3d, cls.mask3d):
            if m:
                f_tot += d * (f - _f) ** 2
                r_tot += d * (r - _r) ** 2
                c_tot += d * (c - _c) ** 2

        _sf = f_tot / d_tot
        _sr = r_tot / d_tot
        _sc = c_tot / d_tot

        cls.goldmasked3dvar = matrix.col((_sf, _sr, _sc))

        # f_tot = 0.0
        # r_tot = 0.0
        # c_tot = 0.0

        # for (f, r, c), d, m in zip(cls.points3d, cls.pixels3d, cls.mask3d):
        #  if m:
        #    f_tot += d * (f - _f) ** 2
        #    r_tot += d * (r - _r) ** 2
        #    c_tot += d * (c - _c) ** 2

        # _sf = f_tot / (d_tot-1)
        # _sr = r_tot / (d_tot-1)
        # _sc = c_tot / (d_tot-1)

        # cls.goldmasked3dubvar = matrix.col((_sf, _sr, _sc))

        pixel_x = flex.double()
        pixel_y = flex.double()
        pixel_z = flex.double()
        pixel_d = flex.double()
        for (x, y, z), d, m in zip(cls.points3d, cls.pixels3d, cls.mask3d):
            if m:
                pixel_x.append(x)
                pixel_y.append(y)
                pixel_z.append(z)
                pixel_d.append(d)

        xc = flex.mean_and_variance(flex.double(pixel_x), pixel_d)
        yc = flex.mean_and_variance(flex.double(pixel_y), pixel_d)
        zc = flex.mean_and_variance(flex.double(pixel_z), pixel_d)
        cls.goldmasked3dubvar = matrix.col(
            (
                xc.gsl_stats_wvariance(),
                yc.gsl_stats_wvariance(),
                zc.gsl_stats_wvariance(),
            )
        )

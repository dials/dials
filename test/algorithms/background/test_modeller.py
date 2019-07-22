from __future__ import absolute_import, division, print_function

import random


class TestExact(object):
    def setup_class(self):
        self.size = (9, 10, 11)

    def test_constant2d_modeller(self):
        from dials.algorithms.background.simple import Constant2dModeller

        modeller = Constant2dModeller()
        eps = 1e-7
        for i in range(10):
            c, data, mask = self.generate_constant_background_2d(self.size, 0, 100)
            model = modeller.create(data, mask)
            assert len(model.params()) == 9
            for j in range(9):
                assert abs(model.params()[j] - c[j]) < eps

    def test_constant3d_modeller(self):
        from dials.algorithms.background.simple import Constant3dModeller

        modeller = Constant3dModeller()
        eps = 1e-7
        for i in range(10):
            c, data, mask = self.generate_constant_background_3d(self.size, 0, 100)
            model = modeller.create(data, mask)
            assert len(model.params()) == 1
            for j in range(1):
                assert abs(model.params()[j] - c) < eps

    def test_linear2d_modeller(self):
        from dials.algorithms.background.simple import Linear2dModeller

        modeller = Linear2dModeller()
        eps = 1e-7
        for i in range(10):
            p, data, mask = self.generate_linear_background_2d(self.size, 0, 100)
            model = modeller.create(data, mask)
            assert len(model.params()) == 3 * 9
            for j in range(9):
                for k in range(3):
                    assert abs(model.params()[k + j * 3] - p[j][k]) < eps

    def test_linear3d_modeller(self):
        from dials.algorithms.background.simple import Linear3dModeller

        modeller = Linear3dModeller()
        eps = 1e-7
        for i in range(10):
            p, data, mask = self.generate_linear_background_3d(self.size, 0, 100)
            model = modeller.create(data, mask)
            assert len(model.params()) == 4
            for j in range(4):
                assert abs(model.params()[j] - p[j]) < eps

    def generate_constant_background_2d(self, size, bmin, bmax):
        from scitbx.array_family import flex

        data = flex.double(flex.grid(size), 0)
        mask = flex.bool(flex.grid(size), True)
        slice_size = (1, size[1], size[2])
        cs = []
        for k in range(size[0]):
            c = random.uniform(bmin, bmax)
            data[k : k + 1, :, :] = flex.double(flex.grid(slice_size), c)
            cs.append(c)
        return cs, data, mask

    def generate_constant_background_3d(self, size, bmin, bmax):
        from scitbx.array_family import flex

        c = random.uniform(bmin, bmax)
        data = flex.double(flex.grid(size), c)
        mask = flex.bool(flex.grid(size), True)
        return c, data, mask

    def generate_linear_background_2d(self, size, bmin, bmax):
        from scitbx.array_family import flex
        from scitbx import matrix

        data = flex.double(flex.grid(size), 0)
        mask = flex.bool(flex.grid(size), True)
        params = []
        for k in range(size[0]):
            a00 = random.uniform(bmin, bmax)
            a01 = random.uniform(bmin, bmax)
            a10 = random.uniform(bmin, bmax)
            p00 = matrix.col((0.5, 0.5, a00))
            p01 = matrix.col((8.5, 0.5, a01))
            p10 = matrix.col((0.5, 8.5, a10))
            n = (p01 - p00).cross(p10 - p00)
            b = n[0]
            c = n[1]
            d = n[2]
            a = -(b * 0.5 + c * 0.5 + d * a00)
            a /= -d
            b /= -d
            c /= -d
            for j in range(size[1]):
                for i in range(size[2]):
                    data[k, j, i] = a + b * (i + 0.5) + c * (j + 0.5)
            eps = 1e-7
            assert abs(data[k, 0, 0] - a00) < eps
            assert abs(data[k, 8, 0] - a10) < eps
            assert abs(data[k, 0, 8] - a01) < eps
            params.append((a, b, c))
        return params, data, mask

    def generate_linear_background_3d(self, size, bmin, bmax):
        from scitbx.array_family import flex
        from scitbx import matrix

        data = flex.double(flex.grid(size), 0)
        mask = flex.bool(flex.grid(size), True)
        a000 = random.uniform(bmin, bmax)
        a001 = random.uniform(bmin, bmax)
        a010 = random.uniform(bmin, bmax)
        a100 = random.uniform(bmin, bmax)
        p000 = matrix.col((0.5, 0.5, 0.5, a000))
        p001 = matrix.col((8.5, 0.5, 0.5, a001))
        p010 = matrix.col((0.5, 8.5, 0.5, a010))
        p100 = matrix.col((0.5, 0.5, 8.5, a100))
        v1 = p001 - p000
        v2 = p010 - p000
        v3 = p100 - p000
        m1 = matrix.sqr((v1[1], v1[2], v1[3], v2[1], v2[2], v2[3], v3[1], v3[2], v3[3]))
        m2 = matrix.sqr((v1[0], v1[2], v1[3], v2[0], v2[2], v2[3], v3[0], v3[2], v3[3]))
        m3 = matrix.sqr((v1[0], v1[1], v1[3], v2[0], v2[1], v2[3], v3[0], v3[1], v3[3]))
        m4 = matrix.sqr((v1[0], v1[1], v1[2], v2[0], v2[1], v2[2], v3[0], v3[1], v3[2]))
        n = matrix.col(
            (m1.determinant(), -m2.determinant(), m3.determinant(), -m4.determinant())
        )
        b = n[0]
        c = n[1]
        d = n[2]
        e = n[3]
        a = -(b * 0.5 + c * 0.5 + d * 0.5 + e * a000)
        a /= -e
        b /= -e
        c /= -e
        d /= -e
        for k in range(size[0]):
            for j in range(size[1]):
                for i in range(size[2]):
                    data[k, j, i] = a + b * (i + 0.5) + c * (j + 0.5) + d * (k + 0.5)
        eps = 1e-7
        assert abs(data[0, 0, 0] - a000) < eps
        assert abs(data[8, 0, 0] - a100) < eps
        assert abs(data[0, 8, 0] - a010) < eps
        assert abs(data[0, 0, 8] - a001) < eps
        return (a, b, c, d), data, mask


class TestPoisson(object):
    def setup_class(self):
        self.size = (9, 9, 9)

    def test_constant2d_modeller(self):
        from dials.algorithms.background.simple import Constant2dModeller
        from dials.array_family import flex

        modeller = Constant2dModeller()
        ma = 10
        sboxes, masks = self.generate_background(self.size, 1000, ma, 0, 0, 0)
        a = []
        v = []
        for i in range(1000):
            model = modeller.create(sboxes[i], masks[i])
            assert len(model.params()) == 9
            assert len(model.variances()) == 9
            a.extend(list(model.params()))
            v.extend(list(model.variances()))

        # Compute Z for each parameter
        z = (flex.double(a) - ma) / flex.sqrt(flex.double(v))

        # Check it looks standard normal
        self.assert_std_norm(z)

    def test_constant3d_modeller(self):
        from dials.algorithms.background.simple import Constant3dModeller
        from dials.array_family import flex

        modeller = Constant3dModeller()

        ma = 10
        sboxes, masks = self.generate_background(self.size, 1000, ma, 0, 0, 0)
        a = []
        v = []
        for i in range(1000):
            model = modeller.create(sboxes[i], masks[i])
            assert len(model.params()) == 1
            assert len(model.variances()) == 1
            a.append(model.params()[0])
            v.append(model.variances()[0])

        # Compute Z for each parameter
        z = (flex.double(a) - ma) / flex.sqrt(flex.double(v))

        # Check it looks standard normal
        self.assert_std_norm(z)

    def test_linear2d_modeller(self):
        from dials.algorithms.background.simple import Linear2dModeller
        from dials.array_family import flex

        modeller = Linear2dModeller()

        # Generate shoeboxes
        ma = 100
        mb = 1
        mc = 2
        sboxes, masks = self.generate_background(self.size, 1000, ma, mb, mc, 0)

        pa = []
        pv = []
        for i in range(1000):
            model = modeller.create(sboxes[i], masks[i])
            assert len(model.params()) == 9 * 3
            assert len(model.variances()) == 9 * 3
            p = model.params()
            v = model.variances()
            for j in range(9):
                pa.append(tuple(p[3 * j : 3 * (j + 1)]))
                pv.append(tuple(v[3 * j : 3 * (j + 1)]))
        a, b, c = zip(*pa)
        va, vb, vc = zip(*pv)

        # Compute Z for each parameter
        za = (flex.double(a) - ma) / flex.sqrt(flex.double(va))
        zb = (flex.double(b) - mb) / flex.sqrt(flex.double(vb))
        zc = (flex.double(c) - mc) / flex.sqrt(flex.double(vc))

        # Check it looks standard normal
        self.assert_std_norm(za)
        self.assert_std_norm(zb)
        self.assert_std_norm(zc)

    def test_linear3d_modeller(self):
        from dials.algorithms.background.simple import Linear3dModeller
        from dials.array_family import flex

        modeller = Linear3dModeller()

        # Generate shoeboxes
        ma = 100
        mb = 1
        mc = 2
        md = 3
        sboxes, masks = self.generate_background(self.size, 1000, ma, mb, mc, md)

        # Compute model
        pa = []
        pv = []
        for i in range(1000):
            model = modeller.create(sboxes[i], masks[i])
            assert len(model.params()) == 4
            assert len(model.variances()) == 4
            p = model.params()
            v = model.variances()
            pa.append(tuple(p))
            pv.append(tuple(v))
        a, b, c, d = zip(*pa)
        va, vb, vc, vd = zip(*pv)

        # Compute Z for each parameter
        za = (flex.double(a) - ma) / flex.sqrt(flex.double(va))
        zb = (flex.double(b) - mb) / flex.sqrt(flex.double(vb))
        zc = (flex.double(c) - mc) / flex.sqrt(flex.double(vc))
        zd = (flex.double(d) - md) / flex.sqrt(flex.double(vd))

        # Check it looks standard normal
        self.assert_std_norm(za)
        self.assert_std_norm(zb)
        self.assert_std_norm(zc)
        self.assert_std_norm(zd)

    def assert_std_norm(self, z):
        from dials.array_family import flex

        mv = flex.mean_and_variance(z)
        m = mv.mean()
        s = mv.unweighted_sample_standard_deviation()
        # Mean-test failure rate:
        # 5*P(abs(X/1000) > 0.1336) where X is normally distributed
        #   with mean 0 and standard deviation sqrt(1000)
        if (abs(m) > 0.1336) or (abs(s - 1.0) > 0.1):
            # from matplotlib import pylab
            # pylab.hist(list(z), 100)
            # pylab.show()
            raise Exception(
                "Mean %f (abs. value <0.1336), Sdev %f (>0.9, <1.1)" % (m, s)
            )

    def generate_background(self, size, N, A, B, C, D):
        from dials.algorithms.simulation.generate_test_reflections import (
            random_background_plane2,
        )
        from dials.array_family import flex
        from dials.util.command_line import ProgressBar

        sboxes = []
        masks = []
        progress = ProgressBar(title="Generating Background")
        for i in range(N):
            mask = flex.bool(flex.grid(size), True)
            sbox = flex.double(flex.grid(size), 0)
            random_background_plane2(sbox, A, B, C, D)
            sboxes.append(sbox)
            masks.append(mask)
            progress.update(100.0 * i / N)
        progress.finished("Generated Background")
        return sboxes, masks

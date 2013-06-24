
class Test(object):

    def __init__(self):
        pass

    def run(self):
        self.tst_getters()
        self.tst_indexing()
        self.tst_nearest()

    def tst_getters(self):
        from math import sqrt
        from dials.algorithms.integration.profile import XdsCircleSampler
        width = 1000
        height = 1000
        depth = 10
        nz = 2
        sampler = XdsCircleSampler((width, height, depth), nz)
        volume_size = sampler.volume_size()
        image_centre = sampler.image_centre()
        r0 = sampler.r0()
        r1 = sampler.r1()
        r2 = sampler.r2()
        size = len(sampler)

        assert(width == volume_size[0])
        assert(height == volume_size[1])
        assert(depth == volume_size[2])
        assert(width // 2 == image_centre[0])
        assert(height // 2 == image_centre[1])
        assert(r0 == min([width // 2, height // 2]))
        assert(r1 == r0 / 3.0)
        assert(r2 == r1 * sqrt(5.0))
        assert(9 * nz == size)
        print 'OK'

    def tst_indexing(self):
        from dials.algorithms.integration.profile import XdsCircleSampler
        from math import sin, cos, pi
        width = 1000
        height = 1000
        depth = 10
        nz = 2
        sampler = XdsCircleSampler((width, height, depth), nz)
        zstep = depth / nz
        eps = 1e-10

        dt = 2 * pi / 8
        r2 = sampler.r2()

        xp = [width / 2.0] + [width / 2.0 + r2 * sin(i * dt) for i in range(8)]
        yp = [height / 2.0] + [height / 2.0 + r2 * cos(i * dt) for i in range(8)]
        zind = [[k] * 9 for k in range(nz)]
        zind = [i for j in zind for i in j]
        zp = [(z + 0.5) * zstep for z in zind]

        for x0, y0, z0, (x1, y1, z1) in zip(xp, yp, zp, sampler):
            assert(abs(x0 - x1) <= eps)
            assert(abs(y0 - y1) <= eps)
            assert(abs(z0 - z1) <= eps)

        print 'OK'


    def tst_nearest(self):
        from math import sqrt, atan2, pi
        from random import randint
        from dials.algorithms.integration.profile import XdsCircleSampler
        width = 1000
        height = 1000
        depth = 10
        nz = 2
        sampler = XdsCircleSampler((width, height, depth), nz)
        xc, yc = sampler.image_centre()
        r1 = sampler.r1()

        for i in range(1000):
            x = randint(0, 1000)
            y = randint(0, 1000)
            z = randint(0, 10)

            r = sqrt((x - xc)**2 + (y - yc)**2)
            if r < r1:
                index00 = 0
            else:
                t = atan2(y - yc, x - xc) + pi
                ai = int(t * 8 / (2 * pi) + 0.5) % 8
                index00 = ai + 1

            index01 = int((z+0.5) * nz // 10)
            index0 = index00 + index01 * 9
            index1 = sampler.nearest((x, y, z))
            assert(index0 == index1)

        print 'OK'


if __name__ == '__main__':
    test = Test()
    test.run()

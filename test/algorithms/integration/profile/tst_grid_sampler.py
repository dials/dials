from __future__ import division

class Test(object):

    def __init__(self):
        pass

    def run(self):
        self.tst_getters()
        self.tst_indexing()
        self.tst_nearest()
        self.tst_self_consistent()

    def tst_getters(self):
        from dials.algorithms.integration.profile import GridSampler
        width = 1000
        height = 1000
        depth = 10
        nx = 10
        ny = 10
        nz = 2
        sampler = GridSampler((width, height, depth), (nx, ny, nz))
        image_size = sampler.image_size()
        grid_size = sampler.grid_size()
        step_size = sampler.step_size()
        size = len(sampler)

        assert(width == image_size[0])
        assert(height == image_size[1])
        assert(depth == image_size[2])
        assert(nx == grid_size[0])
        assert(ny == grid_size[1])
        assert(nz == grid_size[2])
        assert(step_size[0] == width / nx)
        assert(step_size[1] == height / ny)
        assert(step_size[2] == depth / nz)
        assert(nx * ny * nz == size)
        print 'OK'

    def tst_indexing(self):
        from dials.algorithms.integration.profile import GridSampler
        width = 1000
        height = 1000
        depth = 10
        nx = 10
        ny = 10
        nz = 2
        sampler = GridSampler((width, height, depth), (nx, ny, nz))
        xstep, ystep, zstep = sampler.step_size()
        xind = [[i for i in range(nx)]] * ny * nz
        yind = [[j] * nx for j in range(ny)] * nz
        zind = [[k] * nx * ny for k in range(nz)]
        xind = [i for j in xind for i in j]
        yind = [i for j in yind for i in j]
        zind = [i for j in zind for i in j]

        xp = [(x + 0.5) * xstep for x in xind]
        yp = [(y + 0.5) * ystep for y in yind]
        zp = [(z + 0.5) * zstep for z in zind]

        eps = 1e-10

        for x0, y0, z0, (x1, y1, z1) in zip(xp, yp, zp, sampler):
            assert(abs(x0 - x1) <= eps)
            assert(abs(y0 - y1) <= eps)
            assert(abs(z0 - z1) <= eps)

        print 'OK'


    def tst_nearest(self):
        from random import randint
        from dials.algorithms.integration.profile import GridSampler
        width = 1000
        height = 1000
        depth = 10
        nx = 10
        ny = 10
        nz = 2
        sampler = GridSampler((width, height, depth), (nx, ny, nz))

        for i in range(1000):
            x = randint(0, 1000)
            y = randint(0, 1000)
            z = randint(0, 10)
            i = int((x+0.5) * nx // 1000)
            j = int((y+0.5) * ny // 1000)
            k = int((z+0.5) * nz // 10)
            if i >= nx:
                i = nx - 1
            if j >= ny:
                j = ny - 1
            if k >= nz:
                k = nz - 1
            index0 = i + j * nx + k * nx * ny
            index1 = sampler.nearest((x, y, z))
            assert(index0 == index1)

        print 'OK'

    def tst_self_consistent(self):
        from dials.algorithms.integration.profile import GridSampler
        width = 1000
        height = 1000
        depth = 10
        nx = 10
        ny = 10
        nz = 2
        sampler = GridSampler((width, height, depth), (nx, ny, nz))

        for i in range(len(sampler)):
            coord = sampler[i]
            index = sampler.nearest(coord)
            assert(index == i)

        print 'OK'


if __name__ == '__main__':
    test = Test()
    test.run()

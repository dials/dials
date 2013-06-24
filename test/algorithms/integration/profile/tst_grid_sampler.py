
class Test(object):

    def __init__(self):
        pass

    def run(self):
        self.tst_getters()
        self.tst_indexing()
        self.tst_nearest()

    def tst_getters(self):
        from dials.algorithms.integration.profile import GridSampler
        width = 1000
        height = 1000
        nx = 10
        ny = 10
        sampler = GridSampler((width, height), (nx, ny))
        image_size = sampler.image_size()
        grid_size = sampler.grid_size()
        step_size = sampler.step_size()
        size = len(sampler)

        assert(width == image_size[0])
        assert(height == image_size[1])
        assert(nx == grid_size[0])
        assert(ny == grid_size[1])
        assert(step_size[0] == width / nx)
        assert(step_size[1] == height / ny)
        assert(nx * ny == size)
        print 'OK'

    def tst_indexing(self):
        from dials.algorithms.integration.profile import GridSampler
        width = 1000
        height = 1000
        nx = 10
        ny = 10
        sampler = GridSampler((width, height), (nx, ny))
        xstep, ystep = sampler.step_size()
        xind = [[i for i in range(nx)]] * ny
        yind = [[j] * nx for j in range(ny)]
        xind = [i for j in xind for i in j]
        yind = [i for j in yind for i in j]

        xp = [(x + 0.5) * xstep for x in xind]
        yp = [(y + 0.5) * ystep for y in yind]

        eps = 1e-10

        for x0, y0, (x1, y1) in zip(xp, yp, sampler):
            assert(abs(x0 - x1) <= eps)
            assert(abs(y0 - y1) <= eps)

        print 'OK'


    def tst_nearest(self):
        from random import randint
        from dials.algorithms.integration.profile import GridSampler
        width = 1000
        height = 1000
        nx = 10
        ny = 10
        sampler = GridSampler((width, height), (nx, ny))

        for i in range(1000):
            x = randint(0, 1000)
            y = randint(0, 1000)
            i = int((x+0.5) * nx // 1000)
            j = int((y+0.5) * ny // 1000)
            index0 = i + j * nx
            index1 = sampler.nearest((x, y))
            assert(index0 == index1)

        print 'OK'


if __name__ == '__main__':
    test = Test()
    test.run()

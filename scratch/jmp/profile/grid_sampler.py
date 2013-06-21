

class GridSampler(object):

    def __init__(self, width, height, nx, ny):
        self._width = width
        self._height = height
        self._nx = nx
        self._ny = ny
        self._xsize = self._width / self._nx
        self._ysize = self._height / self._ny

    def width(self):
        return self._width

    def height(self):
        return self._height

    def nearest(self, x, y):
        from math import floor
        ix = int(floor(x / self._xsize))
        iy = int(floor(y / self._ysize))

        ix = max([0, ix])
        ix = min([ix, self._width-1])
        iy = max([0, iy])
        iy = min([iy, self._height-1])

        return ix + iy * self._nx

    def __getitem__(self, index):
        x = ((index % self._nx) + 0.5) * self._xsize
        y = ((index // self._nx) + 0.5) * self._ysize
        return x, y

    def __len__(self):
        return self._nx * self._ny

    def __iter__(self):
        for i in range(len(self)):
            yield self.__getitem__(i)


class XdsCircleSampler(object):

    def __init__(self, width, height):
        from math import sqrt
        self._width = width
        self._height = height
        self._xc = width / 2
        self._yc = height / 2
        self._n = 9
        self._r0 = min([self._xc, self._yc])
        self._r1 = self._r0 / 3.0
        self._r2 = self._r1 * sqrt(5.0)

    def width(self):
        return self._width

    def height(self):
        return self._height

    def xc(self):
        return self._xc

    def yc(self):
        return self._yc

    def r0(self):
        return self._r0

    def r1(self):
        return self._r1

    def r2(self):
        return self._r2

    def nearest(self, x, y):
        from math import sqrt, atan2, pi, floor

        xmc = x - self._xc
        ymc = y - self._yc
        r = sqrt(xmc**2 + ymc**2)
        t = atan2(ymc, xmc) + pi

        if r < self._r1:
            return 0

        return int(floor(t * (self._n - 1) / (2 * pi) + 0.5)) % (self._n - 1) + 1

    def __getitem__(self, index):
        from math import pi, sin, cos

        if index < 0 or index >= self._n:
            raise IndexError('Index must be between 0 and 9')

        if index == 0:
            return self._xc, self._yc

        theta = (index - 1)* 2 * pi / (self._n - 1)
        x = self._xc + self._r2 * sin(theta)
        y = self._yc + self._r2 * cos(theta)
        return x, y

    def __len__(self):
        return self._n

    def __iter__(self):
        for i in range(len(self)):
            yield self.__getitem__(i)


sampler = XdsCircleSampler(1000, 1000)

print 'Getting x and y'
x = [xx for xx, yy in sampler]
y = [yy for xx, yy in sampler]

print x, y


#print 'Getting index'
#import numpy
#image = numpy.zeros((1000, 1000), dtype=numpy.int32)
#for j in range(1000):
#    for i in range(1000):
#        image[j,i] = sampler.nearest(i, j)

#print 'Plotting'
#def report_pixel(x, y):
#   v = image[y, x]

#   return "x=%f y=%f value=%f" % (x, y, v)

#from matplotlib import pylab
##ax = pylab.gca()
##ax.format_coord = report_pixel
##pylab.imshow(image, interpolation='none')
#pylab.scatter(x, y)
#pylab.plot()
#pylab.show()

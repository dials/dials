
from dials.algorithms.reflexion_basis.transform import rebin_pixels

class Test(object):

    def __init__(self):
        pass

    def run(self):
        self.tst_identical()
        self.tst_known_offset()
        self.tst_known_orientation()

    def tst_identical(self):
        from scitbx.array_family import flex
        from random import uniform

        # Set the size of the grid
        height = 10
        width = 10

        # Create the grid data
        grid = flex.double([uniform(0, 100) for i in range(height * width)])
        grid.reshape(flex.grid(height, width))

        # Create the grid coordinates
        xy = []
        for j in range(height + 1):
            for i in range(width + 1):
                xy.append((i, j))
        gridxy = flex.vec2_double(xy)
        gridxy.reshape(flex.grid(height+1, width+1))

        # Get the output grid
        output = rebin_pixels(grid, gridxy, (height, width))

        # Check that each pixel is identical
        eps = 1e-7
        for j in range(height):
            for i in range(width):
                assert(abs(output[j,i] - grid[j, i]) <= eps)

        # Test passed
        print 'OK'

    def tst_known_offset(self):
        pass

    def tst_known_orientation(self):
        pass


if __name__ == '__main__':
    test = Test()
    test.run()

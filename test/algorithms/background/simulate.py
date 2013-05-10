

def gaussian(grid, A, cent, sig):

    from scitbx.array_family import flex
    from math import exp
    import numpy

    # Make sure everythings the same dimensions
    assert(len(grid) == len(cent))
    assert(len(grid) == len(sig))

    # Create the grid
    data = flex.double(flex.grid(grid)).as_numpy_array()

    # Calculate the gaussian at each grid point
    for index, value in numpy.ndenumerate(data):
        components = 0
        for x0, s, x in zip(cent, sig, index):
            y = ((x - x0)**2) / (2 * s**2)
            components += y
        data[index] = A * exp(-components)

    # Return the gaussian on the grid
    return flex.double(data)

def poisson_noise(grid, A):
    import numpy
    from scitbx.array_family import flex
    return flex.double(numpy.random.poisson(A, grid).astype(numpy.float64))

if __name__ == '__main__':

    from scitbx.array_family import flex

    spot = gaussian((30,30), 20, (15,15), (3,3))
    background = poisson_noise((30,30), 20)
    grid = background + spot

    from dials.algorithms.background import NormalDiscriminator
    from dials.algorithms.background import PoissonDiscriminator
    from dials.algorithms.background import MeanSubtractor

    discriminate_normal = NormalDiscriminator(n_sigma = 3.0)
    discriminate_poisson = PoissonDiscriminator(n_sigma = 3.0)

    grid_i = flex.int(flex.grid(grid.all()))
    for i, g in enumerate(grid):
        grid_i[i] = int(g)

    mask_normal = discriminate_normal(grid_i)
    mask_poisson = discriminate_poisson(grid_i)

    line = grid.as_numpy_array()[15,:]

    from matplotlib import pylab, cm
    pylab.subplot2grid((2, 3), (0, 0), colspan=3)
    pylab.plot(line)
    pylab.subplot2grid((2, 3), (1, 0))
    pylab.imshow(grid.as_numpy_array(), interpolation='none')
    pylab.subplot2grid((2, 3), (1, 1))
    pylab.imshow(mask_normal.as_numpy_array(), interpolation='none')
    pylab.subplot2grid((2, 3), (1, 2))
    pylab.imshow(mask_poisson.as_numpy_array(), interpolation='none')
    pylab.show()

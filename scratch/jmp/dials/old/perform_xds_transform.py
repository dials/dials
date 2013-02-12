

def extract_reflection_profiles():

    from scitbx import matrix

    import numpy

    from matplotlib import pylab, cm
    from mpl_toolkits.mplot3d import axes3d

    from dials_jmp.old.beam import Beam
    from dials_jmp.old.detector import Detector
    from dials_jmp.old.goniometer import Goniometer
    from dials_jmp.old.reflection.mask import ReflectionMask
    from dials_jmp.old.reflection.grid import Grid
    from dials_jmp.old.reflection.grid_mapper import GridMapper
    from dials_jmp.old.reflection.coordinate_system import CoordinateSystem as ReflectionCoordinateSystem
    from dials_jmp.old.pycbf_ext import search_for_image_volume
    # Set a load of parameters from the GXPARM file
    z0 = 1
    phi0 = 1.0
    dphi = 1.0
    s0 = matrix.col((0.013142, 0.002200, 1.450476))
    dx0 = 244.836136
    dy0 = 320.338531
    px = 0.172000
    py = 0.172000
    d1 = matrix.col((0.000000, -0.939693, -0.342020))
    d2 = matrix.col((1.000000,  0.000000,  0.000000))
    d3 = matrix.col((0.000000, -0.342020,  0.939693))
    m2 = matrix.col((0.999975, -0.001289, -0.006968))
    f  = 122.124901
    l  = 0.689400
    sigma_d = 0.034
    sigma_m = 0.082

    # Set some values of a test reflection
    s1 = matrix.col((-0.0175199028171, -0.24775259748, 1.42844630118))
    phi_dash = 5.0
    phi = 5.83575672475
    sx = 117.588714455
    sy = 311.621428845
    sz = z0 + (phi - phi0) / dphi

    # Load the CBF image volume
    image_volume = search_for_image_volume('/home/upc86896/Projects/data/300k/ximg2700*.cbf')

    image_volume_copy = numpy.copy(image_volume)

    # Make sure vectors are really of length 1/lambda
    s0 = s0.normalize() / l
    s1 = s1.normalize() / l

    # Create the beam, goniometer and detector
    beam = Beam(s0, l)
    gonio = Goniometer(m2, s0, z0, phi0, dphi)
    detector = Detector(image_volume, d1, d2, d3, f, (dx0, dy0), (px, py))

    # Create the reflection coordinate system
    rcs = ReflectionCoordinateSystem(s0, s1, m2, phi)

    rmask = ReflectionMask(sigma_d, sigma_m)

    # Create the grid and mapper to convert detector to reflection coords
    #grid = Grid(1, grid_size=(9, 9, 9), step_size=(0.08, 0.08, 1.0))
    grid = Grid(1, grid_size=(9, 9, 9), sigma_divergence=sigma_d, sigma_mosaicity=sigma_m)
    mapper = GridMapper(gonio, detector, beam, grid)
    mapper.map_reflection(0, rcs, rmask, sx, sy, sz)

#    print 1/0

    from matplotlib import cm, rcParams

    #rcParams.update({'font.size': 6})


    #fig = pylab.figure(figsize=(6,4), dpi=300)
    #for i in range(0,9):
    #    ax = pylab.subplot(3, 3, i+1)
    #    plt = pylab.imshow(image_volume_copy[sz-4-1+i,sy-4:sy+4+1,sx-4:sx+4+1], vmin=0, vmax=2000, cmap=cm.Greys_r)
    #    plt.axes.get_xaxis().set_ticks([])
    #    plt.axes.get_yaxis().set_ticks([])
    #    ax.set_title("slice: {0}".format(i))
    #fig.savefig('/home/upc86896/Documents/Reflection.tiff')
    #pylab.show()


    #fig = pylab.figure(figsize=(6,4), dpi=300)
    #for i in range(0,9):
    #    ax = pylab.subplot(3, 3, i+1)
    #    plt=pylab.imshow(grid.get_grid(0)[i,:,:], vmin=0,
    #        vmax=numpy.max(grid.get_grid(0)), cmap=cm.Greys_r)
    #    plt.axes.get_xaxis().set_ticks([])
    #    plt.axes.get_yaxis().set_ticks([])
    #    ax.set_title("slice: {0}".format(i))
    #fig.savefig('/home/upc86896/Documents/Transformed.tiff')
    #pylab.show()


    # Display the grid
    minx = numpy.min(grid.get_grid_coordinates()[0])
    maxx = numpy.max(grid.get_grid_coordinates()[0])
    miny = numpy.min(grid.get_grid_coordinates()[1])
    maxy = numpy.max(grid.get_grid_coordinates()[1])
    minz = numpy.min(grid.get_grid(0))
    maxz = numpy.max(grid.get_grid(0))
    fig = pylab.figure()
    for i in range(0, 9):
        ax = fig.add_subplot(3, 3, i+1, projection='3d')
        ax.set_xlim([minx, maxx])
        ax.set_ylim([miny, maxy])
        ax.set_zlim([minz, maxz])
        ax.set_autoscalex_on(False)
        ax.set_autoscaley_on(False)
        ax.set_autoscalez_on(False)
        plt=ax.plot_wireframe(grid.get_grid_coordinates()[0][i,:,:],
                          grid.get_grid_coordinates()[1][i,:,:],
                          grid.get_grid(0)[i,:,:])
        ax.set_title("slice {0}".format(i))
#        plt.axes.get_xaxis().set_ticks(map(lambda x: minx + x * (maxx - minx) / 2, range(0, 3)))
#        plt.axes.get_yaxis().set_ticks(map(lambda y: miny + y * (maxy - miny) / 2, range(0, 3)))
#        plt.axes.get_zaxis().set_ticks([])
    #fig.savefig('/home/upc86896/Documents/Wireframe.tiff')
    pylab.show()

#    print grid.grid_data




    # Display the estimated fraction of reflection intensity
    #fig = pylab.figure()
    #for i in range(0, 9):
    #    ax = fig.add_subplot(3, 3, i+1, projection='3d')
    #    ax.plot_wireframe(grid.get_grid_coordinates()[0][i,:,:],
    #                      grid.get_grid_coordinates()[1][i,:,:],
    #                      grid.estimated_reflection_intensity(
    #                            sigma_d, sigma_m)[i,:,:])
    #pylab.show()

if __name__ == '__main__':

    # Extract the reflection profiles
    extract_reflection_profiles()

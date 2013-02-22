from scitbx import matrix
import pycbf
import pycbf_ext
import numpy



if __name__ == '__main__':

    from reflection.grid import Grid
    from reflection.coordinate_system import CoordinateSystem as ReflectionCoordinateSystem
    from detector.coordinate_system import CoordinateSystem as DetectorCoordinateSystem
    from glob import glob

    # Set a load of parameters from the GXPARM file
    z0 = 1
    phi0 = 1.0
    dphi = 1.0
    s0 = matrix.col((0.013142, 0.002200, 1.450476))
    dx0 = 244.836136
    dy0 = 320.338531
    px = 0.172000
    py = 0.172000
    d1 = matrix.col((0.000000, -0.939693, -0.342020)) * px
    d2 = matrix.col((1.000000,  0.000000,  0.000000)) * py
    d3 = matrix.col((0.000000, -0.342020,  0.939693))
    m2 = matrix.col((0.999975, -0.001289, -0.006968))
    f  = 122.124901
    l  = 0.689400
    sigma_d = 0.034
    sigma_m = 0.082

    # Set some values of a test reflection
    s1 = matrix.col((-0.0175198982424, -0.247752596714, 1.42844680608))
    phi_dash = 5.0
    phi = 5.83575672475
    sx = 117.588670459
    sy = 311.621434017
    sz = z0 + (phi - phi0) / dphi
    print "SZ:", sz
#    s1 = matrix.col((0.0530220135471, -0.601982894067, 1.31805682352))
#    phi_dash = 1.0
#    phi = 1.3969438378
#    sx = 301.3044543
#    sy = 346.401674002
#    sz = z0 + (phi - phi0) / dphi

    # Load the CBF image volume
    cbf_path = glob('/home/upc86896/Projects/data/300k/ximg2700*.cbf')
    cbf_path.sort()
    image_volume = pycbf_ext.get_image_volume(cbf_path, 487, 619)

    # Make sure vectors are really of length 1/lambda
    s0 = s0.normalize() / l
    s1 = s1.normalize() / l

    dcs = DetectorCoordinateSystem(d1, d2, d3, f, (dx0, dy0), (px, py))
    rcs = ReflectionCoordinateSystem(s0, s1, m2, phi)


#    grid = Grid(1, grid_size=(9, 9, 9), step_size=(0.03, 0.03, 0.05))
    grid = Grid(1, grid_size=(9, 9, 9), step_size=(0.085, 0.085, 1.0))

    x0 = int(sx) - 4
    x1 = int(sx) + 4
    y0 = int(sy) - 4
    y1 = int(sy) + 4
    z0 = int(sz) - 4
    z1 = int(sz) + 4
    print "Z0:", z0
    print "Z1:", z1

    # Loop through all the detector pixels in the selected imaqes
    CX = numpy.zeros(shape=(y1 - y0 + 1, x1 - x0 + 1), dtype=numpy.float32)
    CY = numpy.zeros(shape=(y1 - y0 + 1, x1 - x0 + 1), dtype=numpy.float32)
    CZ = numpy.zeros(shape=(y1 - y0 + 1, x1 - x0 + 1), dtype=numpy.float32)
    CX_Vals = numpy.zeros(shape=(z1 - z0 + 1, y1 - y0 + 1, x1 - x0 + 1), dtype=numpy.float32)
    CY_Vals = numpy.zeros(shape=(z1 - z0 + 1, y1 - y0 + 1, x1 - x0 + 1), dtype=numpy.float32)
    CZ_Vals = numpy.zeros(shape=(z1 - z0 + 1, y1 - y0 + 1, x1 - x0 + 1), dtype=numpy.float32)
    CI_Vals = numpy.zeros(shape=(z1 - z0 + 1, y1 - y0 + 1, x1 - x0 + 1), dtype=numpy.float32)
    for z in range(z0, z1 + 1):
        for y in range(y0, y1 + 1):
            for x in range(x0, x1 + 1):

                # Get the angle at this z and calculate the coordinate at the
                # given detector coordinate in the profile coordinate frame
                phi_dash = phi0 + (z - 1) * dphi
                s_dash = dcs.to_laboratory_coordinate(x, y)
                s_dash = (1.0 / l) * s_dash.normalize()

                c1, c2, c3 = rcs.from_diffracted_beam_vector(s_dash, phi_dash)

                CX_Vals[z-z0, y-y0, x-x0] = c1
                CY_Vals[z-z0, y-y0, x-x0] = c2
                CZ_Vals[z-z0, y-y0, x-x0] = c3
                CI_Vals[z-z0, y-y0, x-x0] = image_volume[z-1, y, x]

                CX[y-y0, x-x0] = c1
                CY[y-y0, x-x0] = c2
                CZ[y-y0, x-x0] = image_volume[z-1, y, x]

                # Set the values in the profile coordinate grid from the mapped
                # pixel value and distribute values to profile grid
                value = image_volume[z-1, y, x]

                grid.add_point_count(0, (c1, c2, c3), value)



    grid.normalize()


    from matplotlib import pylab, cm
#
#
#    for i in range(0,9):
#        plt = pylab.subplot(3, 3, i+1)
#        plt.imshow(image_volume[int(sz)-1-4+i,sy-4:sy+4+1, sx-4:sx+4+1], vmin=0, vmax=2000, cmap=cm.Greys_r, interpolation='nearest')
#    pylab.show()
#
#    for i in range(0,9):
#        plt = pylab.subplot(3, 3, i+1)
#        plt.imshow(grid[i,:,:], vmin = 0, vmax = numpy.max(grid),cmap=cm.Greys_r, interpolation='nearest')
#
#    pylab.show()
#
#
    from mpl_toolkits.mplot3d import axes3d
#    import matplotlib.pyplot as plt

    fig = pylab.figure()
    for i in range(0, 9):
        ax = fig.add_subplot(3, 3, i+1, projection='3d')

        ax.plot_wireframe(CX_Vals[i,:,:], CY_Vals[i,:,:], CI_Vals[i,:,:])
        ax.set_title('i = {0}; c3 = {1}'.format(i+1, CZ_Vals[i,0,0]))
    pylab.show()


    gx = numpy.arange(9 * 9) % 9
    gy = numpy.arange(9 * 9) / 9
    gx.shape = (9, 9)
    gy.shape = (9, 9)

    fig = pylab.figure()
    for i in range(0, 9):
        ax = fig.add_subplot(3, 3, i+1, projection='3d')

        #ax.plot_wireframe(gx, gy, grid.estimated_reflection_intensity(sigma_d, sigma_m)[i,:,:])#grid.get_grid(0)[i,:,:])
        ax.plot_wireframe(gx, gy, grid.get_grid(0)[i,:,:])
#        ax.plot_wireframe(CX_Vals[i,:,:], CY_Vals[i,:,:], CI_Vals[i,:,:])
        ax.set_title('i = {0}; c3 = {1}'.format(i+1, CZ_Vals[i,0,0]))
    pylab.show()

#    from scipy.interpolate import griddata
#    xpoints = CX_Vals[4,:,:].flat
#    ypoints = CY_Vals[4,:,:].flat
#    values  = CI_Vals[4,:,:].flat
#    points = (xpoints, ypoints)
#
#
##    grid_x = numpy.zeros(shape=(51,51), dtype=numpy.float32)
##    grid_y = numpy.zeros(shape=(51,51), dtype=numpy.float32)
#    grid_x, grid_y = numpy.mgrid[-0.4:0.4:50j, -0.4:0.4:50j]
#
##    for j in range(0, 51+1):
##        for i in range(0, 51+1):
##            grid_x[j,i] = (i - 25) * 0.01
##            grid_y[j,i] = (j - 25) * 0.01
#
#    grid_z2 = griddata(points, values, (grid_x, grid_y), method='cubic')
#
#
#    fig = pylab.figure()
#    ax = fig.add_subplot(1, 1, 1, projection='3d')
#    ax.plot_wireframe(grid_x, grid_y, grid_z2)
##        ax.plot_wireframe(CX_Vals[i,:,:], CY_Vals[i,:,:], CI_Vals[i,:,:,])
##        ax.set_title('i = {0}; c3 = {1}'.format(i+1, CZ_Vals[i,0,0]))
#    pylab.show()

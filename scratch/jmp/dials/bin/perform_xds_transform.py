

def extract_reflection_profiles(paths):
    
    from scitbx import matrix

    import numpy

    from matplotlib import pylab, cm
    from mpl_toolkits.mplot3d import axes3d

    from dials.equipment import Beam
    from dials.equipment import Detector
    from dials.equipment import Goniometer
    from dials.geometry import XdsCoordinateSystem
    from dials.geometry.transform import FromDetectorToBeamVector
    from dials.geometry.transform import FromBeamVectorToXds
    from dials.spot_prediction import ReflectionMask
    from dials.spot_prediction import SubtractBackground
    from dials.geometry.algorithm import XdsTransform
    from dials.geometry.algorithm import XdsTransformGrid
    from dials.io.pycbf_extra import search_for_image_volume
    from scitbx.array_family import flex
    
    # Set a load of parameters from the GXPARM file
    z0 = 1
    phi0 = 1.0
    dphi = 1.0
    s0 = matrix.col((0.013142, 0.002200, 1.450476))
    dx0 = 244.836136
    dy0 = 320.338531
    dx = 487
    dy = 619
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
    image_volume = search_for_image_volume(paths)
    image_volume = flex.int(image_volume)
    image_volume_copy = image_volume.deep_copy()

    # Make sure vectors are really of length 1/lambda
    s0 = s0.normalize() / l
    s1 = s1.normalize() / l

    # Create the beam, goniometer and detector
    beam = Beam(s0, l)
    gonio = Goniometer(m2, phi0, dphi, z0)
    detector = Detector(d1, d2, d3, (dx0, dy0), (px, py), (dx, dy), f)
    
    rcs = XdsCoordinateSystem(s0, s1, m2, phi)
    #dcs = DetectorCoordinateSystem(d1, d2, d3)

    roi_size = (1, 4, 4)
    x0 = int(sx) - roi_size[2]
    x1 = int(sx) + roi_size[2]
    y0 = int(sy) - roi_size[1]
    y1 = int(sy) + roi_size[1]
    z0 = int(sz) - int(gonio.starting_frame) - roi_size[1]
    z1 = int(sz) - int(gonio.starting_frame) + roi_size[1]    

    #reflection_mask = ReflectionMask(image_volume.all(), roi_size=roi_size)
    #reflection_mask.create(flex.vec3_double(1, (sx, sy, sz-z0)))
    
    subtract_background = SubtractBackground(image_volume, roi_size, 0.1, 0.1)
    subtract_background.subtract(flex.vec3_double(1, (sx, sy, sz-gonio.starting_frame)))
    
    xds_grid = XdsTransformGrid(1, (4, 4, 4), sigma_d, sigma_m)
    xds_trans = XdsTransform(xds_grid,
                             image_volume,
                             image_volume.all()[::-1],
                             detector,
                             beam,
                             gonio,
                             roi_size)

    xds_trans.calculate(0, (sx, sy, sz), s1, phi)
    grid = xds_grid.data

    from matplotlib import cm, rcParams
    
#    mask = reflection_mask.mask.as_numpy_array()
    #rcParams.update({'font.size': 6})


    #fig = pylab.figure(figsize=(6,4), dpi=300)
    #for i in range(0,9):
    #    ax = pylab.subplot(3, 3, i+1)
    #    plt = pylab.imshow(image_volume[sz-4-1+i,sy-4:sy+4+1,sx-4:sx+4+1], vmin=0, vmax=2000, cmap=cm.Greys_r)
    #    plt.axes.get_xaxis().set_ticks([])
    #    plt.axes.get_yaxis().set_ticks([])   
    #    ax.set_title("slice: {0}".format(i))     
    #fig.savefig('/home/upc86896/Documents/Reflection.tiff')
    #pylab.show()

    
    #fig = pylab.figure(figsize=(6,4), dpi=300)
    sub_grid = grid[0:1,:,:,:]
    sub_grid.reshape(flex.grid(sub_grid.all()[1:4]))

    print (sub_grid.as_numpy_array() * 100).astype(numpy.int32)

    for i in range(0,9):
        ax = pylab.subplot(3, 3, i+1)
        image = sub_grid[i:i+1,:,:]
        image.reshape(flex.grid(image.all()[1:3]))
        print image.all()
        plt=pylab.imshow(image.as_numpy_array(), vmin=0, 
            vmax=flex.max(sub_grid), cmap=cm.Greys_r)#, interpolation='nearest')
        plt.axes.get_xaxis().set_ticks([])
        plt.axes.get_yaxis().set_ticks([])   
        ax.set_title("slice: {0}".format(i))         
    #fig.savefig('/home/upc86896/Documents/Transformed.tiff')
    pylab.show()


    # Display the grid
    #grid_index = 0
    #minx = numpy.min(grid.get_grid_coordinates()[0])
    #maxx = numpy.max(grid.get_grid_coordinates()[0])
    #miny = numpy.min(grid.get_grid_coordinates()[1])
    #maxy = numpy.max(grid.get_grid_coordinates()[1])
    #minz = numpy.min(grid.get_grid(grid_index))
    #maxz = numpy.max(grid.get_grid(grid_index))
    #fig = pylab.figure()
    #for i in range(0, 9):
    #    ax = fig.add_subplot(3, 3, i+1, projection='3d')
    #    ax.set_xlim([minx, maxx])
    #    ax.set_ylim([miny, maxy])
    #    ax.set_zlim([minz, maxz])
    #    ax.set_autoscalex_on(False)
    #    ax.set_autoscaley_on(False)
    #    ax.set_autoscalez_on(False)
    #    plt=ax.plot_wireframe(grid.get_grid_coordinates()[0][i,:,:], 
    #                      grid.get_grid_coordinates()[1][i,:,:],
    #                      grid.get_grid(grid_index)[i,:,:])
    #    ax.set_title("slice {0}".format(i))
#        plt.axes.get_xaxis().set_ticks(map(lambda x: minx + x * (maxx - minx) / 2, range(0, 3)))
#        plt.axes.get_yaxis().set_ticks(map(lambda y: miny + y * (maxy - miny) / 2, range(0, 3)))
#        plt.axes.get_zaxis().set_ticks([])
    #fig.savefig('/home/upc86896/Documents/Wireframe.tiff')
    #pylab.show()

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

    import sys
    # Extract the reflection profiles
    extract_reflection_profiles(sys.argv[1])    

    
    
    

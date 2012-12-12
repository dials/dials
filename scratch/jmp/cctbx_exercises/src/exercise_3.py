from scitbx import matrix
import pycbf
import pycbf_ext
import numpy



def calculate_profile_coordinate_system_axes(s0, s):
    """Calculate the profile coordinate system axes
    
    :param s0: The incident beam vector
    :param s1: The diffracted beam vector
    :returns: A tuple containing the three axes
    
    """
    e1 = s.cross(s0).normalize()
    e2 = s.cross(e1).normalize()
    e3 = (s + s0).normalize()
    return (e1, e2, e3)


def calculate_diffracted_beam_vector_of_point(x, y, x0, y0, d1, d2, d3, f, l):
    """Calculate diffracted beam vector of detector point
    
    Given a point on the detector surface, calculate the corresponding 
    diffracted beam vector, S, that puts the point onto the surface of the
    Ewald sphere using the equation:
    
        vector = (X - X0)d1 + (Y - Y0)d2 + Fd3
    
    :param x: The x detector coordinate
    :param y: The y detector coordinate
    :param x0: The detector x origin
    :param y0: The detector y origin
    :param d1: The detector x axis vector
    :param d2: The detector y axis vector
    :param d3: The detector normal vector
    :param f: The distance to the detector from the crystal
    :param l: The wavelength, lambda
    :returns: The diffracted beam vector, S_dash
    
    """    
    return (1.0 / l) * ((x - x0) * d1 + (y - y0) * d2 + f * d3).normalize()


def map_to_profile_coordinates(s_dash, s, s0, phi_dash, phi, e1, e2, e3, m2):
    """Map a point defined by the sp vector into the profile coordinate frame.
    
    :param s_dash: The diffracted beam vector corresponding to the point
    :param s1: The diffracted beam vector of the reflection
    :param phi_dash: The real rotation angle
    :param phi1: The rotation angle that puts s1 on the Ewald sphere
    :param e1, e2, e3: The profile coordinate system axes
    :param m2: The rotation axis
    :returns: A tuple containing the mapped coordinates of the sp vector
    
    """    
    from math import pi
    x = s - s0
#    x1 = e1.dot(s_dash - s) / s.length()  
#    x2 = e2.dot(s_dash - s) / s.length()  
#    x3 = e3.dot(x.rotate(m2, phi_dash - phi, True) - x) / x.length()
    x1 = e1.dot(s_dash - s) * (180.0 / (pi * s.length()))  
    x2 = e2.dot(s_dash - s) * (180.0 / (pi * s.length()))  
    x3 = e3.dot(x.rotate(m2, phi_dash - phi, True) - x) * (180.0 / (pi * x.length()))
#    x3 = m2.dot(e1) * (phi_dash - phi)
    return (x1, x2, x3)
   

def detector_coordinate_to_profile_coordinate(x, y, x0, y0, d1, d2, d3, f, l, s, 
                                              s0, phi_dash, phi, e1, e2, e3, m2):
    """Convert a pixel coordinate to a coordinate in the profile frame.
    
    First calculate the diffracted beam vector at the pixel coordinate. Then
    map the diffracted beam vector to the profile coordinates and return the
    components of the e1, e2, e3 exes.
    
    :param x: The x detector coordinate
    :param y: The y detector coordinate
    :param x0: The detector x origin
    :param y0: The detector y origin
    :param d1: The detector x axis vector
    :param d2: The detector y axis vector
    :param d3: The detector normal vector
    :param f: The distance to the detector from the crystal
    :param l: The wavelength, lambda
    :param s1: The diffracted beam vector of the reflection
    :param phi_dash: The real rotation angle
    :param phi1: The rotation angle that puts s1 on the Ewald sphere
    :param e1, e2, e3: The profile coordinate system axes
    :param m2: The rotation axis
    :returns: A tuple containing the mapped coordinates 
    
    """    
    s_dash = calculate_diffracted_beam_vector_of_point(x, y, x0, y0, d1, d2, d3, f, l)
    return map_to_profile_coordinates(s_dash, s, s0, phi_dash, phi, e1, e2, e3, m2)


def calculate_point_weights(w, wi, ws):
    """Calculate the weight to give each grid point.
    
    If the point w is not within the eight grid points, this function will
    give misleading answers.
    
    :param w: The point coordinate
    :param wi: The 8 surrounding grid coordinates
    :param ws: The step size in e1, e2, e3
    :returns: The weight at each surrounding grid point
    
    """
    return map(lambda wj: abs(1.0 - abs(wj[0] - w[0]) / ws[0]) * 
                          abs(1.0 - abs(wj[1] - w[1]) / ws[1]) * 
                          abs(1.0 - abs(wj[2] - w[2]) / ws[2]), wi)


def distribute_pixel_value_to_grid_points(grid, ws, wc, x1, x2, x3, value):
    """Given a pixel value, map the value to the grid. 
    
    This is done by giving a proportion of the pixel value to each of the 
    points surrounding the coordinate (x1, x2, x3) on the grid. Grid spacing is
    assumed to be 1.
    
    :param grid: The grid array
    :param x1: The e1 coordinate
    :param x2: The e2 coordinate
    :param x3: The e3 coordinate
    :param value: The pixel value
    
    """    
    from math import floor
    
    #print "w:", tuple((x1, x2, x3))
    #print "wc:", tuple(wc)
    #print "ws:", tuple(ws)
    
    # The staring indices
    i0 = wc[0] + x1 / ws[0]
    j0 = wc[1] + x2 / ws[1]
    k0 = wc[2] + x3 / ws[2]
    #print "i0, j0, k0", tuple((i0, j0, k0))
    
    if (0 <= i0 < 9-1 and 0 <= j0 < 9-1 and 0 <= k0 < 9-1):

        i0 = floor(i0)
        j0 = floor(j0)
        k0 = floor(k0)    
        i1 = i0 + 1
        j1 = j0 + 1
        k1 = k0 + 1

        #print k0

        #print "I:", i0, "->", i1
        #print "J:", j0, "->", j1
        #print "K:", k0, "->", k1
        
        wi = []
        ind = []
        for k in (k0, k1):
            for j in (j0, j1):
                for i in (i0, i1):
                    wi.append(matrix.col((ws[0]*(i-wc[0]), ws[1]*(j-wc[1]), ws[2]*(k-wc[2]))))
                    ind.append((i, j, k))
    
    
        w = matrix.col((x1, x2, x3))
        pi = calculate_point_weights(w, wi, ws) 
        
        for (i, j, k), p in zip(ind, pi):
            grid[k, j, i] += value * p
    
   
def fraction_of_reflection_intensity(c1, c2, c3, de1, de2, de3, sigma_d, sigma_m):
    from math import exp, sqrt, pi
    
    w1 = de1 * exp((-c1**2) / (2 * sigma_d**2)) / (sqrt(2 * pi) * sigma_d)
    w2 = de2 * exp((-c2**2) / (2 * sigma_d**2)) / (sqrt(2 * pi) * sigma_d)
    w3 = de3 * exp((-c3**2) / (2 * sigma_m**2)) / (sqrt(2 * pi) * sigma_m)
    return w1 * w2 * w3




if __name__ == '__main__':

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
    #image_at_z = image_volume[int(sz)-1,sy-4:sy+4, sx-4:sx+4]
    #print image_volume[int(sz)-1,sy-4:sy+4, sx-4:sx+4]

    # Make sure vectors are really of length 1/lambda
    s0 = s0.normalize() / l
    s1 = s1.normalize() / l

    # Create the profile coordinate system
    e1, e2, e3 = calculate_profile_coordinate_system_axes(s0, s1)
    #print "Profile Coordinate System"
    #print "e1:", tuple(e1), "e1.S0:", e1.dot(s0), "e1.s1:", e1.dot(s1)
    #print "e2:", tuple(e2), "e2.S1:", e2.dot(s1), "e2.e1:", e2.dot(e1)
    #print "e3:", tuple(e3), "e3.e1:", e3.dot(e1), "e3.(s1-s0):", e3.dot(s1-s0)
    
    # The profile grid
    # Set the grid step sizes to 0.1 and origin to the centre
    

    grid = numpy.zeros(shape=(9, 9, 9), dtype=numpy.float32)
#    ws = matrix.col((0.02, 0.02, 0.05))    
    ws = matrix.col((0.1, 0.1, 1.0))   
    wc = matrix.col((int(grid.shape[0]/2), int(grid.shape[1]/2), int(grid.shape[2]/2)))


    # Set the start and end, x, y, z values
    x0 = int(sx) - 2
    x1 = int(sx) + 2
    y0 = int(sy) - 2
    y1 = int(sy) + 2    
    z0 = int(sz) - 2
    z1 = int(sz) + 2
    #z0 = int(sz)
    #z1 = z0
        
    # Loop through all the detector pixels in the selected imaqes
    for z in range(z0, z1 + 1):
        for y in range(y0, y1 + 1):
            for x in range(x0, x1 + 1):
                
                # Get the angle at this z and calculate the coordinate at the
                # given detector coordinate in the profile coordinate frame
                phi_dash = phi0 + (z - 1) * dphi
                c1, c2, c3 = detector_coordinate_to_profile_coordinate(x, y, dx0, dy0, 
                    d1, d2, d3, f, l, s1, s0, phi_dash, phi, e1, e2, e3, m2)

                if x == int(sx) and y == int(sy):
                    print z, c1, c2, c3, phi, phi_dash, image_volume[z-1, y, x]
                #print z, c3, wc[2] + c3 / ws[2]

                # Set the values in the profile coordinate grid from the mapped 
                # pixel value and distribute values to profile grid
                value = image_volume[z-1, y, x]
                distribute_pixel_value_to_grid_points(grid, ws, wc, c1, c2, c3, value)



#    for k in range(0, 9):
#        for j in range(0, 9):
#            for i in range(0, 9):
#                
#                de1 = ws[0]
#                de2 = ws[1]
#                de3 = ws[2]
#                c1 = (i - wc[0]) * ws[0]
#                c2 = (j - wc[1]) * ws[1]
#                c3 = (k - wc[2]) * ws[2]
#                grid[k,j,i] = fraction_of_reflection_intensity(c1, c2, c3, de1, de2, de3, sigma_d, sigma_m)
            
        

#    print image_volume[int(sz)-1-4:int(sz)-1+4+1,sy-4:sy+4+1, sx-4:sx+4+1]
    grid = 100 * grid / numpy.max(grid)
#    print numpy.max(grid)
#    print numpy.array(grid*10000/numpy.max(grid), dtype=numpy.int32)
    print numpy.array(grid, dtype=numpy.int32)


    from matplotlib import pylab, cm

    for i in range(0,9):
        plt = pylab.subplot(3, 3, i+1)
        plt.imshow(image_volume[int(sz)-1-4+i,sy-4:sy+4+1, sx-4:sx+4+1], vmin=0, vmax=2000, cmap=cm.Greys_r, interpolation='nearest')
    pylab.show()

    for i in range(0,9):
        plt = pylab.subplot(3, 3, i+1)
        plt.imshow(grid[i,:,:], vmin = 0, vmax = numpy.max(grid),cmap=cm.Greys_r, interpolation='nearest')
    
    pylab.show()

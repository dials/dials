from cctbx.sgtbx import space_group, space_group_symbols
from scitbx import matrix
from cctbx import uctbx
import pycbf
import numpy




def generate_reflection_indices(uc, dmin):
    """Generate the possible reflection indices from a unit cell object: N.B.
    that these are not symmetry reduced."""
   
    # Get the maximum hkl indices from the unit cell parameters
    maxh, maxk, maxl = uc.max_miller_indices(dmin)

    # The list of reflection indices
    indices = []

    # Loop through the maximum range of hkl indices
    for h in range(-maxh, maxh + 1):
        for k in range(-maxk, maxk + 1):
            for l in range(-maxl, maxl + 1):
                
                # Do not include (0, 0, 0)
                if h == 0 and k == 0 and l == 0:
                    continue
                
                # Ensure the resolution is less than the minimum resolution
                if uc.d((h, k, l)) < dmin:
                    continue
                
                # Append indices to list
                indices.append((h, k, l))

    # Return the list of indices
    return indices


def remove_absent_indices(indices, cell_space_group):
    """From the given list of indices, remove those reflections which should be
    systematic absences according to the given space group."""
   
    # The list of present reflections
    present = []
    
    # Check that each reflection is present
    for hkl in indices:
        if not cell_space_group.is_sys_absent(hkl):
            present.append(hkl)

    # Return those reflections that are present
    return present


def generate_intersection_angles(a_matrix, dmin, wavelength, indices):
    """From an A matrix following the Mosglm convention and the list of indices,
    return a list of phi (hkl) where (typically) there will be 2 records
    corresponding to each h, k, l."""
    
    from rstbx.diffraction import rotation_angles
    from scitbx import matrix
    import math
    
    # Construct an object to calculate the rotation of a reflection about
    # the (0, 1, 0) axis.
    ra = rotation_angles(dmin, a_matrix, wavelength, matrix.col((1, 0, 0)))
        
    # Convert radians to degrees
    r2d = 180.0 / math.pi

    # Loop through all the given reflection indices    
    phi_hkl = []
    for hkl in indices:
        if ra(hkl):
            
            # Calculate the 2 intersection angles for the reflection index
            # and append them to the list of phi-angles
            phis = ra.get_intersection_angles()
            phi_hkl.append((phis[0] * r2d % 360, hkl))
            phi_hkl.append((phis[1] * r2d % 360, hkl))
            
    # Return the list of phi-angles
    return phi_hkl


def select_reflections(phi0, phi1, phi_hkl):
    """Select reflections in range phi0 to phi1 inclusive."""
    return [ph[1] for ph in phi_hkl if (ph[0] >= phi0 and ph[0] <= phi1)]



def generate_observed_reflection_indices(ub_matrix, unit_cell, cell_space_group, 
                                         dmin, wavelength):
    """Predict the reflections.
    
    Calculate the unit cell from the orthogonalization matrix. Then calculate
    the indices of all predicted reflections based on the unit cell parameters
    and the given resolution. Then calculate the symmetry group and remove the
    absent indices. Lastly, get the wavelength from the CBF file and generate
    the intersection angles, phi and return these.
    
    :param cbf_handle: The CBF file handle
    :param ub_matrix: The UB matrix
    :param dmin: The resolution
    :returns: A list of reflection indices
    
    """
    # Generate reflection indices from the unit cell parameters and resolution.
    # Then remove indices that are systemmatically absent because of symmetry.
    # Generate the intersection angles of the remaining reflections and return.
    indices = generate_reflection_indices(unit_cell, dmin)
    present = remove_absent_indices(indices, cell_space_group)
    return generate_intersection_angles(ub_matrix, dmin, wavelength, present)
    
    
def calculate_reflection_detector_coordinates(indices, s0, x0, y0, f, d1, d2, d3, rub):
    """Calculate the pixel coordinates of each reflection.
    
    Calculate the diffracted beam wave vector and find where it intersects
    with the dectector plane.
    
    :param indices: The miller indices
    :param s0: The incident beam wave vector
    :param m2: The crystal rotation axis
    :param b1: The unit cell a axis unit vector
    :param b2: The unit cell b axis unit vector
    :param b3: The unit cell c axis unit vector
    :param x0: The x coordinate of the detector origin
    :param y0: The y coordinate of the detector origin
    :param f: The distance from the crystal to the detector
    :param d1: The detector x axis vector
    :param d2: The detector y axis vector
    :param d3: The detector normal unit vector
    :param phi: The rotation angle
    :returns: A list of pixel indices
    
    """
    # Loop through all the miller indices. Calculate the pixel coordinates of
    # the reflection and append them to the coordinate list.
    coords = []
    for h, k, l in indices:
 
        # Calculate the reciprocal lattice vector and rotate it by phi about
        # the rotation axis. Then calculate the diffracted beam vector, s1.
        x = rub * matrix.col((h, k, l))
        s1 = s0 + x

        # Calculate the point (in pixels) at which the diffracted beam vector 
        # intersects the detector plane. (only when s1.d3 > 0)
        s1_dot_d3 = s1.dot(d3)
        if (s1_dot_d3 > 0):
            x = x0 + f * s1.dot(d1) / s1_dot_d3
            y = y0 + f * s1.dot(d2) / s1_dot_d3
            coords.append((x, y))

    # Return the pixel coordinates
    return coords


def read_image_from_cbf(cbf_handle, category='array_data', column='data', 
                        row=0, element=0, as_3d=False):
    """Read an image from a CBF file
    
    This function is a bit of a hack - I'm not sure what the general structure
    of a CBF file is like but for the data I have, it works. Reads an image
    from the location specified in the CBF file, otherwise raises an exception.
    
    :param cbf_handle: The handle to the CBF file
    :param category: Category in which the image is contained
    :param column: Column in which the image is contained
    :param row: Row in which image is contained
    :param element: Element in which image is contained
    :param as_3d: Reshape the array in either 2D or 3D
    :returns: An array of image data
    
    """
    # Find the given category, column and row
    cbf_handle.find_category(category)
    cbf_handle.find_column(column)
    cbf_handle.select_row(row)
    
    # Check the type of the element to ensure it's a binary
    # otherwise raise an exception
    type = cbf_handle.get_typeofvalue()
    if type.find('bnry') > -1:

        # Read the image data into an array
        image_string = cbf_handle.get_integerarray_as_string()
        image = numpy.fromstring(image_string, numpy.int32)

        # Get the size of the image data (either 2d or 3d)
        if as_3d == True:
            image_size = cbf_handle.get_3d_image_size(element)
        else:
            image_size = cbf_handle.get_image_size(element)
        
        # Resize the image
        image.shape = (image_size)

    else:
        raise TypeError('{0}:{1}:{2}:{3} is not an image'.format(
                            category, column, row, element))

    # Return the image
    return image


def read_reflections_from_cbf(cbf_handle, positions, bbox = (5,5)):
    """Read the reflections from the CBF file.
    
    For each reflection specified in the positions array - which should be of
    the form [(x0, y0), (x1, y1), ...] - extract an portion of the image
    like (x0-bbox[0]:x0+bbox[0], y0-bbox[1]:y0+bbox[1]). Return an array of
    arrays containing the extracted reflection images.
    
    :param cbf_handle: The handle to the cbf file
    :param positions: The array of reflection positions
    :param bbox: The box coordinates for the reflection sub-images
    :returns: An array of reflection images of the size specified by bbox
    
    """
    from matplotlib import pylab, cm
   
    # Read the image from the cbf file
    image = read_image_from_cbf(cbf_handle)

    # Get the width and height of the image
    height, width = image.shape

    pylab.imshow(image, vmin=0, vmax=1000, cmap = cm.Greys_r, origin='lower')
    
    good_coords = [coord for coord in positions if coord[0] >= 0 and coord[1] >= 0 and coord[0] < width and coord[1] < height]
    x = [i for i, j in good_coords]
    y = [j for i, j in good_coords]
    #x = [i for i, j in positions]
    #y = [j for i, j in positions]
    print "Len:", len(x), len(y)
    print "Min, Max Image:", numpy.min(image), numpy.max(image)
    print "Min X, Y:", numpy.min(x), numpy.min(y)    
    print "Max X, y:", numpy.max(x), numpy.max(y)
    
    pylab.scatter(x, y, marker='+')
    pylab.show()
    
    # bounding box distances.
    bx, by = bbox
    
    # For each reflection position, extract a region around the
    # reflection of the size given by bbox, append this sub-image
    # into the reflection profile array
    profiles = []
    for (i, j) in positions:
        
        # The sub-image coordinates
        i0, j0, i1, j1 = (i - bx, j - by, i + bx + 1, j + by + 1)
        
        # Ensure the position (including the bounding box) is within the image, 
        # if not, continue to the next reflection position
        if (i0 <= 0 or j0 <= 0 or i1 >= width or j1 >= height):
            continue

        # Get the sub-image
        profiles.append(image[i0:i1, j0:j1])

    # Return profile images
    return profiles

def extract_reflections(cbf_path, hdf_path, ub_matrix, bbox, dmin):
    """Extract the reflections from a CBF file and save to HDF5.

    From the supplied orientation matrix do the following:
     - predict all the reflections that will occur
     - compute the positions on the reflections on the detector
     - select a 3D volume, defined by the bbox parameter, from the image 
       around each reflection in the CBF file
     - write these reflection profiles to a hdf5 file    
    
    :param cbf_path: The path to the CBF file
    :param hdf_path: The path to the HDF file
    :param ub_matrix: The UB orientation matrix
    :param bbox: The profile box parameters
    :param dmin: The resolution
    
    """
    
    # Create the cbf file handle and read the file
    cbf_handle = pycbf.cbf_handle_struct()
    cbf_handle.read_file(cbf_path, pycbf.MSG_DIGEST)
    
    phi0 = 1.0
    dphi = 1.0
    phi1 = phi0 + dphi
    
    wavelength = 0.689400
    pixel_size = (0.172000, 0.172000)
    m2 = matrix.col((0.999975, -0.001289, -0.006968))
    s0 = matrix.col((0.013142, 0.002200, 1.450476))
    f = 122.124901
    r = m2.axis_and_angle_as_r3_rotation_matrix(1.0, deg=True)
    rub = r * ub_matrix
    
    d1 = matrix.col((0.000000, -0.939693, -0.342020)) / pixel_size[0]
    d2 = matrix.col((1.000000,  0.000000,  0.000000)) / pixel_size[1]
    d3 = matrix.col((0.000000, -0.342020,  0.939693))
    dx0, dy0 = 244.836136, 320.338531
    
    # Calculate the unit cell and space group
    unit_cell = uctbx.unit_cell(parameters = (51.5784, 51.5784, 51.5784, 90.000, 90.000, 90.000))
    #unit_cell = uctbx.unit_cell(orthogonalization_matrix=ub_matrix.inverse())
       
    #cell_space_group = unit_cell.lattice_symmetry_group()
    cell_space_group = space_group(space_group_symbols(211).hall())
     
    # Generate the reflections
    reflections = generate_observed_reflection_indices(ub_matrix, unit_cell, 
        cell_space_group, dmin, wavelength)

    # Select the reflections on this frame
    reflections = select_reflections(phi0, phi1, reflections)
    
    # Calculate the reflection detector coordinates    
    coords = calculate_reflection_detector_coordinates(
        reflections, s0, dx0, dy0, f, d1, d2, d3, rub)
    
    # Read the reflection profiles from the CBF file
    profiles = read_reflections_from_cbf(cbf_handle, coords, bbox)

    # Write the reflection profiles to a HDF5 file    
    #write_reflections_to_hdf5(hdf_path, profiles, cbf_path, ub_matrix, dmin)
    

if __name__ == '__main__':
    
    # The CBF image path
    cbf_path = '../data/ximg2700_00001.cbf'
    
    # Set the HDF file path
    hdf_path = '../data/ximg2701_00001_reflection_profiles.hdf5'
    
    #A = UB
    #ub_matrix = matrix.sqr([ 0.0144873, -0.0128813, -0.0002988,
    #                        -0.0128113, -0.0143530, -0.0024004,
    #                         0.0013736,  0.0019910, -0.0192366])
    
    iub_matrix = matrix.sqr([33.915398, -36.200676, -14.127665,
                             38.855934,  31.328993,  13.001726,
                             -0.544143, -19.192179,  47.871685])
 
    ub_matrix = iub_matrix.inverse()
    
    # Set the size of the reflection profile box
    bbox = (5, 5)
    
    # Set the resolution
    dmin = 0.7
   
    # From the supplied orientation matrix do the following:
    #  - predict all the reflections that will occur
    #  - compute the positions on the reflections on the detector
    #  - select a 3D volume from the image around each reflection (in cbf file)
    #  - write these reflection profiles to a hdf5 file
    extract_reflections(cbf_path, hdf_path, ub_matrix, bbox, dmin)
    

#    
#    from xdsio import IntegrateFile
#    integrate_handle = IntegrateFile()
#    integrate_handle.read_file('../data/INTEGRATE.HKL')
#    xyzobs             = integrate_handle.xyzobs
#    xyzcal             = integrate_handle.xyzcal
#
#    xcoord = []
#    ycoord = []
#    for x, y, z in xyzcal:
#        if (z >= 1 and z < 2):
#            xcoord.append(x)
#            ycoord.append(y)
#    
#    cbf_handle = pycbf.cbf_handle_struct()
#    cbf_handle.read_file(cbf_path, pycbf.MSG_DIGEST)    
#    image = read_image_from_cbf(cbf_handle)
#    height, width = image.shape
#
#    from matplotlib import pylab, cm
#    pylab.imshow(image, vmin=0, vmax=1000, cmap = cm.Greys_r, origin='lower')
#    
#    print numpy.max(xcoord), numpy.max(ycoord)
#    
#    pylab.scatter(xcoord, ycoord, marker='+')
#    pylab.show()    
    
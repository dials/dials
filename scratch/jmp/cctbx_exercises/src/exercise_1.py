import logging
import os
import numpy
import pycbf
import h5py
from scitbx import matrix
from cctbx import uctbx
from cctbx.sgtbx import space_group


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
    logging.debug('Finding image in cbf file...')

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


def generate_reflection_indices(uc, dmin):
    """Generate the possible reflection indices from a unit cell object: N.B.
    that these are not symmetry reduced."""
    logging.debug('Generating reflection indices...')

    # Get the maximum hkl indices from the unit cell parameters
    maxh, maxk, maxl = uc.max_miller_indices(dmin)

    print maxh, maxk, maxl

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
    logging.debug('Removing absent indices...')

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


def generate_observed_reflection_indices(cbf_handle, ub_matrix, dmin):
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
    logging.debug('Generating Reflections...')

    # Calculate the unit cell parameters from the orthogonalization matrix
    unit_cell = uctbx.unit_cell(orthogonalization_matrix=ub_matrix.inverse())
    logging.debug('Calculated unit cell parameters:\n'
                  '\ta={0}\n\tb={1}\n\tc={2}\n'
                  '\talpha={3}\n\tbeta={4}\n\tgamma={5}'.format(
                        *(unit_cell.parameters())))

    # Generate reflection indices
    indices = generate_reflection_indices(unit_cell, dmin)
    logging.debug('Generated {0} reflection indices'.format(len(indices)))

    # Calculate the symmetry group
    cell_space_group = unit_cell.lattice_symmetry_group()
    logging.debug('Calculated space group: {0}'.format(cell_space_group.info()))

    # Remove absent indices
    present = remove_absent_indices(indices, cell_space_group)
    logging.debug('Removed {0} absent reflections'.format(
                    len(indices)-len(present)))

    # Get the wavelength from the file
    wavelength = cbf_handle.get_wavelength()
    logging.debug('Got wavelength from file: {0} A'.format(wavelength))

    # Generate the intersection angles
    observable = generate_intersection_angles(ub_matrix, dmin, wavelength, present)
    logging.debug('Generated {0} intersection angles.'.format(len(observable)))

    # Get the phi angles from the CBF file
    dphi = 1.0
    phi0 = 1.0
    phi1 = phi0 + dphi

    # Select the observable reflections within the given phi-angles
    present = select_reflections(phi0, phi1, observable)

    # Return the intersection angles
    return present



def find_beam_direction(cbf_handle):
    """Find the beam direction (why is this not simpler in pycbf?)"""
    cbf_handle.find_category('axis')
    cbf_handle.find_column('equipment')
    cbf_handle.find_row('source')
    beam_direction = []
    for j in range(3):
        cbf_handle.find_column('vector[%d]' % (j + 1))
        beam_direction.append(cbf_handle.get_doublevalue())

    B = - matrix.col(beam_direction).normalize()
    return B


def compute_central_rotation_matrix(gonio):
    """Compute the central rotation matrix"""
    x = gonio.rotate_vector(0.5, 1, 0, 0)
    y = gonio.rotate_vector(0.5, 0, 1, 0)
    z = gonio.rotate_vector(0.5, 0, 0, 1)
    R = matrix.sqr(x + y + z).transpose()
    return R


def calculate_spot_detector_positions(ub_matrix, wavelength, beam_direction,
                                      detector_distance, detector_normal,
                                      rotation_matrix, miller_indices):
    """Calculate the position on the detector of the reflections

    The reflections are specified by their miller_indices. In order to find the
    detector coordinate of the spot, the geometry of the Ewald sphere is used.

    The reciprocal space vector to the reflection can be found by using the
    UB and R matrices using the equation x = RUBh. The vector, S, can be found
    by S = S0 + x, where S0 is the incident beam vector.

    The position on the detector can then be found by finding the intersection
    of the line defined by S, starting at the origin of the crystal with the
    detector plane.

    The detector plane is specified by its normal and a distance from the
    detector to the crystal. The detector origin can then be found as
    dorig = distance * normal

    :param ub_matrix: The UB matrix
    :param wavelength: The wavelength of the incident beam
    :param beam_direction: The beam direction
    :param detector_distance: The distance from the detector to the crystal
    :param detector_normal: The normal of the detector plane
    :param rotation_matrix: The rotation matrix
    :param miller_indices: The miller indices of the reflections
    :returns: The spot detector coordinates

    """
    # Setup parameters for equations below
    RUB = rotation_matrix * ub_matrix
    S0 = (1.0 / wavelength) * beam_direction
    F = detector_distance
    N = matrix.col(detector_normal)

    # Calculate the detector origin
    Dorig = F * (N)
    DorigN = Dorig.dot(N)

    # For each hkl in the set of miller indices calculate the vector x followed
    # by the reciprocal space vector to the point, S. The detector position can
    # then be found using simple geometry to find the point where the vector, S,
    # intersections the detector plane (defined by the detector normal and the
    # detector origin.
    detector_positions = []
    for h, k, l in miller_indices[0:10]:
        X = RUB * matrix.col((h, k, l))
        S = S0 + X
        print S0.length(), S.length()
        SN = S.dot(N)
        if (SN != 0):
            D = DorigN / SN
            P = D * S
            detector_positions.append(P.elems)

    # Return the detector positions
    return detector_positions


def calculate_spot_pixel_indices(cbf_handle, detector, detector_positions):
    """Calculate the pixel coordinates of the detector positions and filter
    to ensure only good indices within the limits of the detector are used.

    :param cbf_handle: The handle to the CBF file
    :param detector: The handle to the detector in the CBF file
    :param detector_positions: The detector positions of the reflections
    :returns: The pixel indices.

    """
    # Get the width and height of the detector
    width, height = cbf_handle.get_image_size_fs(0)

    # Get the coordinate of the pixel at index (0, 0). Then get the
    # coordinates of the pixels at index (0, height-1) and (width-1, 0) and
    # find the vector between them and the origin (at (0, 0)).
    O = matrix.col(detector.get_pixel_coordinates_fs(0, 0))
    DX = matrix.col(detector.get_pixel_coordinates_fs(width-1,  0)) - O
    DY = matrix.col(detector.get_pixel_coordinates_fs(0, height-1)) - O

    # Calculate the length of A and B
    len_DX = DX.length()
    len_DY = DY.length()

    # For each detector position, get the vector between the position and the
    # detector (0, 0) pixel and calculate the length of the vector. Find the
    # pixel coordinate, (i, j) for the reflection by calculating the distance
    # along the detector x and y, DX and DY, vectors.
    pixel_indices = []
    for x, y, z in detector_positions:
        DP = matrix.col((x ,y, z)) - O
        len_DP = DP.length()
        if (len_DP > 0):
            costheta1 = DX.dot(DP) / (len_DX * len_DP)
            costheta2 = DY.dot(DP) / (len_DY * len_DP)
            i = len_DP * costheta1 * (width-1)  / len_DX
            j = len_DP * costheta2 * (height-1) / len_DY
        else:
            i = 0
            j = 0
        pixel_indices.append((i, j))

    # Check each pixel index and make sure that they are within the image limits
    good_pixel_indices = []
    for i, j in pixel_indices:
        if i >= 0 and j >= 0 and i < width and j < height:
            good_pixel_indices.append((i,j))

    # Return the good pixel indices
    return good_pixel_indices


def calculate_reflection_pixel_indices(cbf_handle, reflections, ub_matrix):
    """Calculate the reflection pixel indices on the detector.

    First Calculate the x, y, z coordinate of the reflection in lab coordinates
    where it intersects with the detector plane. Then calculate the pixel
    indices on the detector. Filter the indices so that only valid indices
    are returned.

    :param cbf_handle: The handle to the CBF file
    :param reflections: The list of reflections
    :param ub_matrix: The UB matrix
    :returns: The list of reflection pixel indices

    """
    logging.debug('Computing reflection positions...')

    # Construct the detector object from the CBF file
    detector = cbf_handle.construct_detector(0)

    # Calculate the detector normal, distance and pixel size
    detector_normal = tuple(detector.get_detector_normal())
    detector_distance = detector.get_detector_distance()
    logging.debug("Detector normal: {0}".format(detector_normal))
    logging.debug("Detector distance: {0}".format(detector_distance))

    # Find the beam direction
    beam_direction = find_beam_direction(cbf_handle)
    logging.debug('Beam Direction: {0}'.format(beam_direction))

    # Construct the goniometer and get the rotation matrix
    gonio = cbf_handle.construct_goniometer()
    rotation_matrix = compute_central_rotation_matrix(gonio)
    logging.debug('Rotation Matrix: {0}'.format(rotation_matrix))

    # Get the wavelength
    wavelength = cbf_handle.get_wavelength()
    logging.debug('Got wavelength from file: {0} A'.format(wavelength))

    # Calculate the position on the detector of the reflections
    detector_positions = calculate_spot_detector_positions(ub_matrix,
        wavelength, beam_direction, detector_distance, detector_normal,
        rotation_matrix, reflections)
    logging.debug('Calculated {0} detector positions...'.format(
                    len(detector_positions)))

    # Calculate the pixel indices of the reflections
    pixel_indices = calculate_spot_pixel_indices(cbf_handle, detector,
                                                 detector_positions)
    logging.debug('Calculated {0} pixel indices...'.format(
                    len(pixel_indices)))

    # Return the pixel indices
    return pixel_indices


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
    logging.debug('Reading reflection profiles from CBF file...')

    from matplotlib import pylab, cm

    # Read the image from the cbf file
    image = read_image_from_cbf(cbf_handle)

    pylab.imshow(image, vmin=0, vmax=1000, cmap = cm.Greys_r)

    x = [i for i, j in positions]
    y = [j for i, j in positions]
    print image[173,309], image[315,7],image[203,212]
    print "Min, Max Image:", numpy.min(image), numpy.max(image)
    print "Min X, Y:", numpy.min(x), numpy.min(y)
    print "Max X, y:", numpy.max(x), numpy.max(y)

    pylab.scatter(x, y, marker='+')
    pylab.show()


    # Get the width and height of the image
    width, height = image.shape

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


def write_reflections_to_hdf5(hdf_path, profiles, cbf_path, ub_matrix, dmin):
    """Save all the reflection profiles in a HDF5 file

    :param hdf_path: The path to the HDF5 file
    :param profiles: The reflection profiles
    :param cbf_path: The path to the cbf path from which the data came
    :param ub_matrix: The UB matrix used to calculate the reflections
    :param dmin: The resolution used to calculate the reflections

    """
    logging.debug('Writing reflection profiles to HDF5 file...')

    # Open the HDF5 file
    h5_handle = h5py.File(hdf_path, 'w')

    # Create a header group and a number of datasets containing the header info
    group = h5_handle.create_group('header')
    group.create_dataset('cbf_file', data=os.path.basename(cbf_path))
    group.create_dataset('ub_matrix', data=ub_matrix)
    group.create_dataset('dmin', data=dmin)

    # Create a group to contain the reflection profiles. Then loop through
    # the reflection profiles and save each in a seperate dataset.
    group = h5_handle.create_group('reflection_profiles')
    for idx, val in enumerate(profiles):
        group.create_dataset('{0}'.format(idx), data=val)

    # Close the hdf file
    h5_handle.close()


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
    logging.debug('Opening CBF file...')
    cbf_handle = pycbf.cbf_handle_struct()
    cbf_handle.read_file(cbf_path, pycbf.MSG_DIGEST)

    # Generate the reflections
    reflections = generate_observed_reflection_indices(cbf_handle, ub_matrix, dmin)

    # Compute the reflection positions
    indices = calculate_reflection_pixel_indices(cbf_handle, reflections, ub_matrix)

    # Read the reflection profiles from the CBF file
    profiles = read_reflections_from_cbf(cbf_handle, indices, bbox)

    # Write the reflection profiles to a HDF5 file
    write_reflections_to_hdf5(hdf_path, profiles, cbf_path, ub_matrix, dmin)



if __name__ == '__main__':

    # Set the logging level to show debug messages
    logging.basicConfig(level=logging.DEBUG)

    # The CBF image path
    cbf_path = '../data/ximg2701_00001.cbf'

    # Set the HDF file path
    hdf_path = '../data/ximg2701_0001_reflection_profiles.hdf5'

    # A = UB
    #ub_matrix = matrix.sqr([ 0.0144873, -0.0128813, -0.0002988,
    #                        -0.0128113, -0.0143530, -0.0024004,
    #                         0.0013736,  0.0019910, -0.0192366])

#    print ub_matrix.inverse()

    ub_matrix = matrix.sqr([ 0.0127485,  0.0146061, -0.0002041,
                            -0.0136081,  0.0117765, -0.0072142,
                            -0.0053110,  0.0048870,  0.0179949])

    # Set the size of the reflection profile box
    bbox = (5, 5)

    # Set the resolution
    dmin = 0.7

    # Print some information about the CBF file
#    import cbfinfo
#    cbfinfo.print_info(cbf_path)

    # From the supplied orientation matrix do the following:
    #  - predict all the reflections that will occur
    #  - compute the positions on the reflections on the detector
    #  - select a 3D volume from the image around each reflection (in cbf file)
    #  - write these reflection profiles to a hdf5 file
    extract_reflections(cbf_path, hdf_path, ub_matrix, bbox, dmin)

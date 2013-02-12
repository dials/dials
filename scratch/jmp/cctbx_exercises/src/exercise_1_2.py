"""This module is a test exercise using CCTBX to do the spot prediction in the
XDS method and to save the resulting spot profiles to a HDF file."""

from cctbx.sgtbx import space_group, space_group_symbols
from scitbx import matrix
from cctbx import uctbx
import pycbf
import h5py
import xdsio
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


def generate_observed_reflection_indices(ub_matrix, unit_cell, cell_space_group,
                                         dmin, wavelength):
    """Predict the reflections.

    Calculate the indices of all predicted reflections based on the unit cell
    parameters and the given resolution. Then remove the indices that are
    absent due to the symmetry. Lastly, generate the intersection angles, phi,
    and return these.

    :param cbf_handle: The CBF file handle
    :param ub_matrix: The UB matrix
    :param unit_cell: The unit cell parameters
    :param cell_space_group: The unit cell space group symmetry
    :param dmin: The resolution
    :param wavelength: The beam wavelength
    :returns: A list of reflection indices

    """
    # Generate reflection indices from the unit cell parameters and resolution.
    # Then remove indices that are systemmatically absent because of symmetry.
    # Generate the intersection angles of the remaining reflections and return.
    indices = generate_reflection_indices(unit_cell, dmin)
    present = remove_absent_indices(indices, cell_space_group)
    return generate_intersection_angles(ub_matrix, dmin, wavelength, present)


def calculate_reflection_detector_coordinates(indices, angles, m2, s0, x0, y0,
                                              f, d1, d2, d3, ub):
    """Calculate the pixel coordinates of each reflection.

    Calculate the diffracted beam wave vector and find where it intersects
    with the dectector plane.

    :param indices: The miller indices
    :param angles: The phi angles of the reflections
    :param m2: The rotation axis
    :param s0: The incident beam wave vector
    :param x0: The x coordinate of the detector origin
    :param y0: The y coordinate of the detector origin
    :param f: The distance from the crystal to the detector
    :param d1: The detector x axis vector
    :param d2: The detector y axis vector
    :param d3: The detector normal unit vector
    :param ub: The UB matrix
    :returns: A list of pixel coordinates

    """
    # Loop through all the miller indices. Calculate the pixel coordinates of
    # the reflection and append them to the coordinate list.
    coords = []
    index1 = 0
    index2 = 0
    for hkl, phi in zip(indices, angles):

        # Calculate the reciprocal lattice vector and rotate it by phi about
        # the rotation axis. Then calculate the diffracted beam vector, s1.
        rub = m2.axis_and_angle_as_r3_rotation_matrix(phi, deg=True) * ub
        x = rub * matrix.col(hkl)
        s1 = s0 + x

        # Calculate the point (in pixels) at which the diffracted beam vector
        # intersects the detector plane. (only when s1.d3 > 0)
        s1_dot_d3 = s1.dot(d3)
        if (s1_dot_d3 > 0):
            x = x0 + f * s1.dot(d1) / s1_dot_d3
            y = y0 + f * s1.dot(d2) / s1_dot_d3
            coords.append((x, y))
            if (int(x) == 117 and int(y) == 311):
                print index2, x, y, s1, phi

            index1 += 1

        index2 += 1

    # Return the pixel coordinates
    return coords


def read_image_from_cbf(cbf_handle, category='array_data', column='data',
                        row=0, element=0):
    """Read an image from a CBF file

    This function is a bit of a hack - I'm not sure what the general structure
    of a CBF file is like but for the data I have, it works. Reads an image
    from the location specified in the CBF file, otherwise raises an exception.

    :param cbf_handle: The handle to the CBF file
    :param category: Category in which the image is contained
    :param column: Column in which the image is contained
    :param row: Row in which image is contained
    :param element: Element in which image is contained
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
        image_size = cbf_handle.get_image_size(element)

        # Resize the image
        image.shape = (image_size)

    else:
        raise TypeError('{0}:{1}:{2}:{3} is not an image'.format(
                            category, column, row, element))

    # Return the image
    return image


def read_reflections_from_volume(volume, coords, bbox = (5,5,5)):
    """Read the reflections from the CBF file.

    For each reflection specified in the positions array - which should be of
    the form [(x0, y0, z0), (x1, y1, z1), ...] - extract an portion of the
    image volume like:
        (x0-bbox[0]:x0+bbox[0], y0-bbox[1]:y0+bbox[1], z0-bbox[2]:z0+bbox[2]).
    Return an array of arrays containing the extracted reflection images.

    :param volume: The image volume
    :param coords: The array of reflection positions
    :param bbox: The box coordinates for the reflection sub-images
    :returns: An array of reflection images of the size specified by bbox

    """
    # Get the width and height of the image
    sz, sy, sx = volume.shape

    # bounding box distances.
    bx, by, bz = bbox

    # For each reflection position, extract a region around the
    # reflection of the size given by bbox, append this sub-image
    # into the reflection profile array
    profiles = []
    for (i, j, k) in coords:

        # The sub-image coordinates
        i0, i1 = i - bx, i + bx + 1
        j0, j1 = j - by, j + by + 1
        k0, k1 = k - bz, k + bz + 1

        # Ensure the position (including the bounding box) is within the image,
        # if not, continue to the next reflection position
        if (i0 < 0 or j0 < 0 or k0 < 0 or i1 >= sx or j1 >= sy or k1 >= sz):
            continue

        # Get the sub-image, numpy uses [col, row] format
        profiles.append(volume[k0:k1, j0:j1, i0:i1])

    # Return profile images
    return profiles


def select_reflections(phi0, phi1, phi_hkl):
    """Select reflections in range phi0 to phi1 inclusive."""
    return [ph for ph in phi_hkl if (ph[0] >= phi0 and ph[0] <= phi1)]


def load_cbf_image_volume(cbf_paths, width, height):
    """Load the image volume from the list of cbf_paths. The list of paths is
    assumed to be is order from 1->n.

    :param cbf_paths: The list of cbf files
    :param width The width (xsize) of the volume
    :param height The height (ysize) of the volume
    :returns: The 3D volume array

    """
    # Initialise the image volume
    num_slices = len(cbf_paths)
    volume = numpy.zeros(shape=(num_slices, height, width), dtype=numpy.int32)

    # For each CBF file, read the image and put into the image volume
    for i, filename in enumerate(cbf_paths):
        cbf_handle = pycbf.cbf_handle_struct()
        cbf_handle.read_file(filename, pycbf.MSG_DIGEST)
        volume[i,:,:] = read_image_from_cbf(cbf_handle)

    # Return the image volume
    return volume


def write_reflections_to_hdf5(hdf_path, volume, coords, profiles):
    """Save all the reflection profiles in a HDF5 file

    :param hdf_path: The path to the HDF5 file
    :param volume: The image volume
    :param coords: The reflection coordinates
    :param profiles: The reflection profiles

    """
    # Open the HDF5 file
    h5_handle = h5py.File(hdf_path, 'w')

    # Create a group to hold image volume
    group = h5_handle.create_group('image_volume')
    group.create_dataset('data', data=volume)

    # Create a group with reflection positions
    group = h5_handle.create_group('reflections')
    group.create_dataset('detector_coordinates', data=coords)

    # Create a group to contain the reflection profiles. Then loop through
    # the reflection profiles and save each in a seperate dataset.
    group = h5_handle.create_group('reflection_profiles')
    for idx, val in enumerate(profiles):
        group.create_dataset('{0}'.format(idx), data=val)

    # Close the hdf file
    h5_handle.close()


def extract_and_save_reflections(cbf_path, gxparm_path, hdf_path, bbox, dmin):
    """Extract the reflections from a CBF file and save to HDF5.

    Do the following:
     - Read all the data from the given input files
     - predict all the reflections that will occur
     - compute the positions on the reflections on the detector
     - select a 3D volume, defined by the bbox parameter, from the image
       around each reflection in the CBF file
     - write these reflection profiles to a hdf5 file

    :param cbf_path: The path to the CBF file
    :param gxparm_path: The path to the GXPARM file
    :param hdf_path: The path to the HDF file
    :param bbox: The profile box parameters
    :param dmin: The resolution

    """
    # Create the GXPARM file and read the contents
    gxparm_handle = xdsio.GxParmFile()
    gxparm_handle.read_file(gxparm_path)

    # Read the wavelength, pixel size, detector origin and distance
    wavelength = gxparm_handle.wavelength
    pixel_size = gxparm_handle.pixel_size
    detector_size = gxparm_handle.detector_size
    dx0, dy0 = gxparm_handle.detector_origin
    f = gxparm_handle.detector_distance

    # Read the space group symmetry and unit cell parameters
    symmetry = gxparm_handle.space_group
    unit_cell = uctbx.unit_cell(parameters = gxparm_handle.unit_cell)
    cell_space_group = space_group(space_group_symbols(symmetry).hall())

    # Read the starting angle, angle delta and starting frame
    phi0 = gxparm_handle.starting_angle
    dphi = gxparm_handle.oscillation_range
    z0   = gxparm_handle.starting_frame

    # Read the rotation axis and beam vector
    m2 = matrix.col(gxparm_handle.rotation_axis)
    s0 = matrix.col(gxparm_handle.beam_vector).normalize() / wavelength

    # Read the unit cell and detector axis vectors
    b1 = matrix.col(gxparm_handle.unit_cell_a_axis)
    b2 = matrix.col(gxparm_handle.unit_cell_b_axis)
    b3 = matrix.col(gxparm_handle.unit_cell_c_axis)
    d1 = matrix.col(gxparm_handle.detector_x_axis) / pixel_size[0]
    d2 = matrix.col(gxparm_handle.detector_y_axis) / pixel_size[1]
    d3 = matrix.col(gxparm_handle.detector_normal)

    # Create the UB matrix
    ub_matrix = matrix.sqr(b1.elems + b2.elems + b3.elems).inverse()

    # Load the image volume from the CBF files
    volume = load_cbf_image_volume(cbf_path, detector_size[0], detector_size[1])
    volume_size_z, volume_size_y, volume_size_x = volume.shape

    # Generate the reflections. Get all the indices at the given resolution.
    # Then using the space group symmetry, remove any indices that will be
    # systemmatically absent. Finally, calculate the intersection angles of
    # each of the reflections. The returned variable contains a list of
    # (phi, hkl) elements
    phi_hkl = generate_observed_reflection_indices(ub_matrix, unit_cell,
        cell_space_group, dmin, wavelength)

    # Calculate the minimum and maximum angles and filter the reflections to
    # leave only those whose intersection angles lie between them
    phi_min = phi0
    phi_max = phi0 + (volume_size_z - z0) * dphi
    phi_hkl = select_reflections(phi_min, phi_max, phi_hkl)
    hkl = [h for p, h in phi_hkl]
    phi = [p for p, h in phi_hkl]

    # Calculate the reflection detector coordinates. Calculate the
    # diffracted beam vector for each reflection and find the pixel
    # coordinate where the line defined by the vector intersects the
    # detector plane. Returns a list of (x, y) detector coordinates.
    coords = calculate_reflection_detector_coordinates(hkl, phi, m2, s0, dx0,
        dy0, f, d1, d2, d3, ub_matrix)

    # Calculate the z coordinate (frame number) of each reflection and create
    # an array containing a list of (x, y, z) coordinates
    zcoords = [(angle - phi0) / dphi + z0 for angle in phi]
    coords = [(x, y, z-1) for (x, y), z in zip(coords, zcoords)]

    index1 = 0
    index2 = 0
    for x, y, z in coords:
        if 0 <= x < volume_size_x and 0 <= y < volume_size_y and 0 <= z < volume_size_z:
            if (index1 == 97663):
                print index2
            index1 += 1
        index2 += 1

    # Filter the coordinates to those within the boundaries of the volume
    coords = [(x, y, z) for x, y, z in coords
                if 0 <= x < volume_size_x and
                   0 <= y < volume_size_y and
                   0 <= z < volume_size_z]



    # Read the reflections from the volume. Return a 3D profile of each
    # reflection with a size of (2*bbox[0]+1, 2*bbox[1]+1, 2*bbox[2]+1)
    profiles = read_reflections_from_volume(volume, coords, bbox)

    # Write the reflection profiles to a HDF5 file
    write_reflections_to_hdf5(hdf_path, volume, coords, profiles)

    # Print points on first image
    from matplotlib import pylab, cm

    image_size = volume.shape
    image_size_ratio = image_size[2] / float(image_size[1])
    figure_size = (6, 6 / image_size_ratio)
    for i in range(0, 1):
        #fig = pylab.figure(figsize=figure_size, dpi=300)
        image = volume[i,:,:]
        filtered_xy = [(x, y) for x, y, z in coords if i <= z < i+1]
        xcoords = [x for x, y in filtered_xy]
        ycoords = [y for x, y in filtered_xy]
        intensities = [image[y, x] for x, y in filtered_xy]
        #index = numpy.where(intensities == numpy.max(intensities))[0][0]
        #mean_intensity = numpy.mean(intensities)
        #diff_mean = numpy.abs(numpy.array(intensities)-mean_intensity)
        #index = numpy.where(diff_mean == numpy.min(diff_mean))[0][0]
    #    print index, xcoords[index], ycoords[index], intensities[index], numpy.mean(intensities)
        plt = pylab.imshow(image, vmin=0, vmax=1000, cmap=cm.Greys_r)
        pylab.scatter(xcoords, ycoords, marker='x')
        plt.axes.get_xaxis().set_ticks([])
        plt.axes.get_yaxis().set_ticks([])
        #fig.savefig('/home/upc86896/Documents/image_volume_w_markers{0}.tiff'.format(i),
        #    bbox_inches='tight', pad_inches=0)
        pylab.show()

if __name__ == '__main__':

    from glob import glob

    # The CBF image path, make sure paths are in order
    cbf_path = glob('/home/upc86896/Projects/data/300k/ximg2700*.cbf')
    cbf_path.sort()

    # Set the GXPARM path
    gxparm_path = '../data/GXPARM.XDS'

    # Set the HDF file path
    hdf_path = '../data/ximg2700_reflection_profiles.hdf5'

    # Set the size of the reflection profile box
    bbox = (5, 5, 5)

    # Set the resolution
    dmin = 0.7

    # From the supplied orientation matrix do the following:
    #  - predict all the reflections that will occur
    #  - compute the positions on the reflections on the detector
    #  - select a 3D volume from the image around each reflection (in cbf file)
    #  - write these reflection profiles to a hdf5 file
    extract_and_save_reflections(cbf_path, gxparm_path, hdf_path, bbox, dmin)

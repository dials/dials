#!/usr/bin/env python

"""This module is a test exercise using CCTBX to do the spot prediction in the
XDS method and to save the resulting spot profiles to a HDF file."""

from cctbx.sgtbx import space_group, space_group_symbols, space_group_type
from scitbx import matrix
from cctbx import uctbx
from dials_jmp.old.lattice_point import LatticePoint
from dials_jmp.old.reflection.reflection import Reflection
from dials_jmp.old.index_generator import IndexGenerator
from dials_jmp.io import pycbf_extra
#import h5py
from dials_jmp.io import xdsio
import numpy
from dials_jmp.equipment import Beam, Goniometer, Detector
from dials_jmp.geometry.transform import FromHklToBeamVector
from dials_jmp.geometry.transform import FromHklToDetector
from dials_jmp.geometry.transform import FromBeamVectorToDetector

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
    beam = gxparm_handle.get_beam()
    gonio = gxparm_handle.get_goniometer()
    detector = gxparm_handle.get_detector()
    ub_matrix = gxparm_handle.get_ub_matrix()

    # Read the space group symmetry and unit cell parameters
    symmetry = gxparm_handle.space_group

    #cell_space_group = space_group(space_group_symbols(symmetry).hall())
    cell_space_group_type = space_group_type(space_group(space_group_symbols(symmetry).hall()))

    # Load the image volume from the CBF files
    volume = pycbf_extra.search_for_image_volume(cbf_path)
    volume_size_z, volume_size_y, volume_size_x = volume.shape

    print "Generate Reflections"

    from dials_jmp.spot_prediction import SpotPredictor
    n_image_frames = volume_size_z
    gonio.num_frames = n_image_frames

    uc = uctbx.unit_cell(orthogonalization_matrix = ub_matrix)

    spot_predictor = SpotPredictor(beam, gonio, detector, ub_matrix, dmin,
                                   uc, cell_space_group_type)

    spot_predictor.predict_spots()
    xyz = spot_predictor.image_volume_coordinates
    print "Num Valid: ", len(xyz)


    coords = xyz

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
        #intensities = [image[y, x] for x, y in filtered_xy]
       #index = numpy.where(intensities == numpy.max(intensities))[0][0]
        #mean_intensity = numpy.mean(intensities)
        #diff_mean = numpy.abs(numpy.array(intensities)-mean_intensity)
        #index = numpy.where(diff_mean == numpy.min(diff_mean))[0][0]
        #print index, xcoords[index], ycoords[index], intensities[index], numpy.mean(intensities)
        plt = pylab.imshow(image, vmin=0, vmax=1000, cmap=cm.Greys_r)
        pylab.scatter(xcoords, ycoords, marker='x')
        plt.axes.get_xaxis().set_ticks([])
        plt.axes.get_yaxis().set_ticks([])
       #fig.savefig('/home/upc86896/Documents/image_volume_w_markers{0}.tiff'.format(i),
        #    bbox_inches='tight', pad_inches=0)
        pylab.show()

def test():

    from glob import glob

    # The CBF image path, make sure paths are in order
    #cbf_path = r'C:\Users\upc86896\Documents\Projects\data\300k\ximg2700*.cbf'
    cbf_path = r'/home/upc86896/Projects/data/300k/ximg2700*.cbf'

    # Set the GXPARM path
    gxparm_path = r'C:\Users\upc86896\Documents\Projects\dials\dials-svn\dials-code\scratch\jmp\cctbx_exercises\data\GXPARM.XDS'
    gxparm_path = r'/home/upc86896/Projects/dials/dials-svn/dials-code/scratch/jmp/data/GXPARM.XDS'


    # Set the HDF file path
    hdf_path = '/home/upc86896/Projects/dials/dials-svn/dials-code/scratch/jmp/data/ximg2700_reflection_profiles.hdf5'

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

if __name__ == '__main__':
    test()
    #import cProfile
    #cProfile.run('test()', 'profile.txt')
    #import pstats
    #p = pstats.Stats('profile.txt')
    #p.sort_stats('cumulative').print_stats(10)

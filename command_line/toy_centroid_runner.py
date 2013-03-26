#
# toy_centroid_runner.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: Graeme Winter
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import division

class ToyCentroidRunner(object):
    '''Class to run the centroid script.'''

    def __init__(self, xparm_file, integrate_file, image_files,
                 output_file=None, algorithm=None, selection=None):
        '''Initialise the script.'''

        # Set the required input arguments
        self.xparm_file = xparm_file
        self.integrate_file = integrate_file
        self.image_files = image_files

        # Set the optional output/algorithm and selection
        self.output_file = output_file
        self.algorithm = algorithm
        self.selection = selection

    def __call__(self):
        '''The main body of the script.'''

        # Begin by loading models from the input files
        self._load_models()

        # Get the miller indices from the integrate.hkl file
        self._get_miller_indices()

        # Predict the reflections from the given miller indices
        self._predict_reflections()

        # Extract reflection profiles
        self._extract_reflection_profiles()

        # Calculate centroids
        self._calculate_centroids()

        # Output data to screen and file
        self._write_output()

    def _load_models(self):
        '''Load the models from file.'''
        from iotbx.xds import xparm, integrate_hkl
        from dxtbx.sweep import SweepFactory
        from dials.util import io
        import dxtbx

        # Load the sweep from the image files
        print "Reading image files for sweep."
        self.sweep = SweepFactory.sweep(self.image_files)

        # Load the models from the xparm file
        print "Reading: \"{0}\"".format(self.xparm_file)
        models = dxtbx.load(self.xparm_file)
        self.beam = models.get_beam()
        self.detector = models.get_detector()
        self.gonio = models.get_goniometer()
        self.scan = models.get_scan()

        # Set the scan model with the image range
        image_range = self.sweep.get_scan().get_image_range()
        self.scan.set_image_range(image_range)

        # Read other data (need to assume an XPARM file)
        xparm_handle = xparm.reader()
        xparm_handle.read_file(self.xparm_file)
        self.UB = io.get_ub_matrix_from_xparm(xparm_handle)
        self.unit_cell = io.get_unit_cell_from_xparm(xparm_handle)
        self.space_group = io.get_space_group_type_from_xparm(xparm_handle)
        self.d_min = self.detector.get_max_resolution_at_corners(
            self.beam.get_direction(), self.beam.get_wavelength())

        # Read the integrate file to get the sigma_d and sigma_m
        print "Reading: \"{0}\"".format(self.integrate_file)
        self.integrate_handle = integrate_hkl.reader()
        self.integrate_handle.read_file(self.integrate_file)
        self.sigma_divergence = self.integrate_handle.sigma_divergence
        self.sigma_mosaicity = self.integrate_handle.sigma_mosaicity

        print ""
        print "Experimental Models"
        print "-------------------"
        print self.beam
        print self.detector
        print self.gonio
        print self.scan
        print self.UB

    def _get_miller_indices(self):
        '''Read the miller indices from the integrate.hkl file.'''
        from cctbx.array_family import flex

        # For each record in the integrate.hkl file get the miller index
        self.miller_indices = flex.miller_index()
        for hkl, I, sigI in zip(self.integrate_handle.hkl,
                                self.integrate_handle.iobs,
                                self.integrate_handle.sigma):

            # Select strong reflections
            if sigI <= 0 or I / sigI < 40:
                continue

            # append the miller indices
            self.miller_indices.append(hkl)

    def _predict_reflections(self):
        '''Predict the reflection positions and bounding boxes.'''
        from dials.algorithms.spot_prediction import IndexGenerator
        from dials.algorithms.spot_prediction import RayPredictor
        from dials.algorithms.spot_prediction import ray_intersection
        from dials.algorithms.spot_prediction import reflection_frames
        from dials.algorithms.integration import BBoxCalculator
        from math import pi

        # Create the spot predictor
        s0 = self.beam.get_s0()
        m2 = self.gonio.get_rotation_axis()
        UB = self.UB
        dphi = self.scan.get_oscillation_range(deg=False)
        predict_rays = RayPredictor(s0, m2, UB, dphi)

        # Predict the reflections
        self.reflections = reflection_frames(self.scan, ray_intersection(
            self.detector, predict_rays(self.miller_indices)))

        # Set the divergence and mosaicity
        n_sigma = 5.0
        delta_divergence = n_sigma * self.sigma_divergence * pi / 180.0
        delta_mosaicity = n_sigma * self.sigma_mosaicity * pi / 180.0

        # Create the bounding box calculator
        calculate_bbox = BBoxCalculator(self.beam, self.detector, self.gonio,
            self.scan, delta_divergence, delta_mosaicity)

        # Calculate the frame numbers of all the reflections
        calculate_bbox(self.reflections)

    def _extract_reflection_profiles(self):
        '''Extract the reflection profiles.'''
        from dials.algorithms.integration import extract_reflection_profiles

        # extract reflection profiles
        extract_reflection_profiles(self.sweep, self.reflections)

    def _calculate_centroids(self):
        '''Calculate the centroids of the reflections.'''
        from dials.algorithms.centroid import toy_centroid
        from dials.algorithms.centroid import toy_centroid_lui
        from dials.algorithms.centroid import FilteredCentroid

        # Create the algorithm corresponding to the given option
        if self.algorithm == 'toy':
            tc = toy_centroid(self.reflections)
        elif self.algorithm == 'lui':
            tc = toy_centroid_lui(self.reflections)
        elif self.algorithm == 'xds':
            tc = FilteredCentroid(self.reflections)
        else:
            raise ValueError('No centroid algorithm named {0}'.format(
                self.algorithm))

        # Get the reflections
        self.reflections = tc.get_reflections()

    def _write_output(self):
        '''Write output to std::out and files.'''
        from reflection_stats import ReflectionStats
        import pickle

        # Select a certain number of random reflections
        self._select_random_reflections()

        # Print the reflections
        #for ref in self.reflections:
        #    print ref

        # Print some reflection statistics
        print ReflectionStats(self.reflections)

        # Dump the reflections to file
        if self.output_file:
            print "\nPickling reflection list."
            pickle.dump(self.reflections, open(self.output_file, 'wb'))

    def _select_random_reflections(self):
        '''Select a number of random reflections.'''
        from dials.model.data import ReflectionList
        from random import sample

        # If a selection has been set
        if self.selection:

            # Get a random sample of num indices
            indices = sample(range(len(self.reflections)), self.selection)

            # Create a new reflection list and put selected reflections in
            reflections_new = ReflectionList()
            for i in indices:
                reflections_new.append(self.reflections[i])

            # Return the new list of reflections
            self.reflections = reflections_new


def algorithm_callback(option, opt, value, parser):
    """Parse centroid algorithm"""
    if value != "toy" and value != "lui" and value != "xds":
        raise ValueError("Invalid centroid algorithm.")
    setattr(parser.values, option.dest, value)

if __name__ == '__main__':

    from optparse import OptionParser

    # Specify the command line options
    usage = "usage: %prog [options] /path/to/GXPARM.XDS"
    usage += " /path/to/INTEGRATE.HKL /path/to/image.cbf"

    # Create an option parser
    parser = OptionParser(usage)

    # Add command line options
    parser.add_option(
        '-o', '--output-file',
        dest = 'output_file', type = "string", default = "",
        help = 'Enter a destination filename for reflections')
    parser.add_option(
        '-a', '--algorithm',
        dest = 'algorithm', type = "string", default = "toy",
        action="callback", callback=algorithm_callback,
        help = 'Enter the centroid algorithm to use (toy, lui, xds)')
    parser.add_option(
        '-n', '--num-reflections',
        dest = 'selection', type = "int", default = 0,
        help = 'Enter the number of reflections to select')

    # Parse the arguments
    options, args = parser.parse_args()

    # Print help if no arguments specified, otherwise call function
    if len(args) < 3:
        print parser.print_help()
    else:
        #toy_centroid_runner(args[0], args[1], args[2:], options.output_file,
        #    options.algorithm, options.n_select)

        # Initialise the centroid runner
        centroid_runner = ToyCentroidRunner(
            xparm_file=args[0],
            integrate_file=args[1],
            image_files=args[2:],
            output_file=options.output_file,
            algorithm=options.algorithm,
            selection=options.selection)

        # Run the script
        centroid_runner()

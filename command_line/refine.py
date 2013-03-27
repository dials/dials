#!/usr/bin/env python
#
# refine.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import division

class RefinementRunner(object):
    '''Class to run the refinement script.'''

    def __init__(self, xparm_file, integrate_file, image_files):
        '''Initialise the script.'''

        # Set the required input arguments
        self.xparm_file = xparm_file
        self.integrate_file = integrate_file
        self.image_files = image_files

    def __call__(self):
        '''The main body of the script.'''

        # Begin by loading models from the input files
        self._load_models()

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

def read_reflection_file(filename):
    '''Read reflections from pickle file.'''
    from dials.model.data import ReflectionList
    import pickle
    return pickle.load(open(reflection_file, 'rb'))

def run(reflection_file):
    '''Do the refinement.'''

    # Load the reflections from the pickle file
    reflections = read_reflection_file(reflection_file)

    # Print some reflections
    for r in reflections:
        print r

if __name__ == '__main__':
    from optparse import OptionParser

    # Specify the command line options
    usage = "usage: %prog [options] /path/to/data.p "
    usage += "/path/to/XPARM.XDS "
    usage += "/path/to/INTEGRATE.HKL "
    usage += "/path/to/image*.cbf"

    # Parse the command line options
    parser = OptionParser(usage)
    options, args = parser.parse_args()

    # Print help if no arguments specified, otherwise call spot prediction
    if len(args) < 4:
        print parser.print_help()
    else:

        # Get stuff from args
        reflection_file = args[0]
        xparm_file = args[1]
        integrate_file = args[2]
        image_files = args[3:]

        # Run the refinement
#        run(reflection_file)
        runner = RefinementRunner(xparm_file, integrate_file, image_files)
        runner()

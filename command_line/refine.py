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

    def __init__(self, xparm_file, integrate_file,
                 reflections):
        '''Initialise the script.'''

        # Set the required input arguments
        self.xparm_file = xparm_file
        self.integrate_file = integrate_file
        self.reflections = reflections

    def __call__(self):
        '''The main body of the script.'''
        from scitbx import matrix
        from math import sqrt

        # Begin by loading models from the input files
        self._load_models()

        # pull out data needed for refinement
        temp = [(ref.miller_index, ref.rotation_angle,
            matrix.col(ref.beam_vector), ref.image_coord_mm,
            ref.centroid_variance) for ref in self.reflections]

        hkls, angles, svecs, intersects, variances = zip(*temp)

        # tease out tuples to separate lists
        d1s, d2s = zip(*intersects)
        var_d1s, var_d2s, var_angles = zip(*variances)

        px_size = self.detector.get_pixel_size()
        im_width = self.scan.get_oscillation(deg=False)[1]

        # warning! hard coded nastiness!
        sweep_range = (0, 900 * im_width)

        # change variances to sigmas and convert units
        sig_d1s = [px_size[0] * sqrt(e) for e in var_d1s]
        sig_d2s = [px_size[1] * sqrt(e) for e in var_d2s]
        sig_angles = [im_width * sqrt(e) for e in var_angles]

        #for i in range(10):
        #    print hkls[i], svecs[i], d1s[i], sig_d1s[i], d2s[i], sig_d2s[i], angles[i], sig_angles[i]
        assert len(hkls) == len(svecs) == len(d1s) == len(d2s) == \
               len(sig_d2s) == len(angles) == len(sig_angles)

        from dials.scratch.dgw.refinement import print_model_geometry, refine
        from dials.scratch.dgw.crystal_model import Crystal

        # build a Crystal
        # need a_vec, b_vec, c_vec in the lab frame

        from rstbx.cftbx.coordinate_frame_converter import \
            coordinate_frame_converter
        cfc = coordinate_frame_converter(self.xparm_file)

        a_vec = cfc.get_c('real_space_a')
        b_vec = cfc.get_c('real_space_b')
        c_vec = cfc.get_c('real_space_c')

        mycrystal = Crystal(a_vec, b_vec, c_vec)

        print "Prior to refinement the experimental model is:"
        print_model_geometry(self.beam, self.detector, mycrystal)

        refine(self.beam, self.gonio, mycrystal, self.detector, im_width,
               sweep_range, hkls, svecs, d1s, sig_d1s, d2s, sig_d2s, angles,
               sig_angles, verbosity = 1)

        print
        print "Refinement has completed with the following geometry:"
        print_model_geometry(mybeam, mydetector, mycrystal)

    def _load_models(self):
        '''Load the models from file.'''
        from iotbx.xds import xparm, integrate_hkl
        from dxtbx.sweep import SweepFactory
        from dials.util import io
        import dxtbx

        # Load the models from the xparm file
        print "Reading: \"{0}\"".format(self.xparm_file)
        models = dxtbx.load(self.xparm_file)
        self.beam = models.get_beam()
        self.detector = models.get_detector()
        self.gonio = models.get_goniometer()
        self.scan = models.get_scan()

        # Set the scan model with the image range
        image_range = (1,900)
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
    from scitbx import matrix

    # Load the reflections from the pickle file
    reflections = read_reflection_file(reflection_file)



if __name__ == '__main__':
    from optparse import OptionParser

    # Specify the command line options
    usage = "usage: %prog [options] /path/to/data.p "
    usage += "/path/to/XPARM.XDS "
    usage += "/path/to/INTEGRATE.HKL "

    # Parse the command line options
    parser = OptionParser(usage)
    options, args = parser.parse_args()

    # Print help if no arguments specified, otherwise call spot prediction
    if len(args) < 3:
        print parser.print_help()
    else:

        # Get stuff from args
        reflection_file = args[0]
        xparm_file = args[1]
        integrate_file = args[2]

        # reconstitute the reflections
        reflections = read_reflection_file(reflection_file)
        # Run the refinement
#        run(reflection_file)
        runner = RefinementRunner(xparm_file, integrate_file,
                                  reflections)
        runner()

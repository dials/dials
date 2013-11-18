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
               reflections, verbosity=0, scan_varying=False):
    '''Initialise the script.'''

    # Set the required input arguments
    self.xparm_file = xparm_file
    self.integrate_file = integrate_file
    self.reflections = reflections
    self.verbosity = verbosity
    self.scan_varying = scan_varying

  def __call__(self):
    '''The main body of the script.'''
    from scitbx import matrix
    from math import sqrt
    import random

    # Begin by loading models from the input files
    self._load_models()

    self._saved_reflections = self.reflections.deep_copy()

    # check that the beam vectors are stored: if not, compute them
    # perhaps we should just be computing them here...

    for ref in self.reflections:
      if ref.beam_vector != (0.0, 0.0, 0.0):
        continue
      x, y = self.detector.millimeter_to_pixel(ref.image_coord_mm)
      ref.beam_vector = matrix.col(self.detector.get_pixel_lab_coord(
          (x, y))).normalize() / self.beam.get_wavelength()

    from dials.algorithms.refinement.refinement_helpers import \
        print_model_geometry, refine
    random.seed(42)
    print "Random seed set to 42\n"

    print "Prior to refinement the experimental model is:"
    print_model_geometry(self.beam, self.detector, self.crystal)

    refine(self.beam, self.gonio, self.crystal, self.detector,
           self.scan, self.reflections, verbosity = self.verbosity, fix_cell=False,
           scan_varying=self.scan_varying)

    print
    print "Refinement has completed with the following geometry:"
    print_model_geometry(self.beam, self.detector, self.crystal)

    # Do a test of new reflection pos
    self._update_reflections_test()

  def _load_models(self):
    '''Load the models from file.'''
    from iotbx.xds import xparm, integrate_hkl
    from dxtbx.imageset import ImageSetFactory
    from dials.util import ioutil
    import dxtbx
    from rstbx.cftbx.coordinate_frame_converter import \
        coordinate_frame_converter
    from scitbx import matrix

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
    print "\nWARNING: the image range has been hardcoded to (1,900)"
    print "If that does not match the input data, this script will"
    print "probably fail!\n"

    # Read other data (need to assume an XPARM file)
    xparm_handle = xparm.reader()
    xparm_handle.read_file(self.xparm_file, check_filename = False)
    self.space_group_type = ioutil.get_space_group_type_from_xparm(xparm_handle)
    cfc = coordinate_frame_converter(self.xparm_file)
    a_vec = cfc.get_c('real_space_a')
    b_vec = cfc.get_c('real_space_b')
    c_vec = cfc.get_c('real_space_c')
    self.unit_cell = cfc.get_unit_cell()
    self.UB = matrix.sqr(a_vec.elems + b_vec.elems + c_vec.elems).inverse()

    from dials.model.experiment.crystal_model import Crystal
    self.crystal = Crystal(a_vec, b_vec, c_vec,
                           space_group = self.space_group_type.group())

    # Calculate resolution
    d_min = self.detector.get_max_resolution(self.beam.get_s0())

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
    print self.crystal

  def _update_reflections_test(self):
    from cctbx.array_family import flex
    from collections import defaultdict

    # Get miller indices from saved reflectons
    miller_indices = [r.miller_index for r in self._saved_reflections]

    self.miller_indices = flex.miller_index(miller_indices)

    print "Predicting new reflections"
    self._predict_reflections()

    # Put coords from same hkl in dict for saved reflections
    coord1 = defaultdict(list)
    for r1 in self._saved_reflections:
      coord1[r1.miller_index].append(r1.image_coord_px)

    # Put coords from same hkl in dict for new reflections
    coord2 = defaultdict(list)
    for r2 in self._new_reflections:
      coord2[r2.miller_index].append(r2.image_coord_px)

    # Print out coords for each hkl
    for h in coord1.keys():
      c1 = coord1[h]
      c2 = coord2[h]
      #print c1, c2

  def _predict_reflections(self):
    '''Predict the reflection positions and bounding boxes.'''
    from dials.algorithms.spot_prediction import IndexGenerator
    from dials.algorithms.spot_prediction import RayPredictor
    from dials.algorithms.spot_prediction import ray_intersection
    from dials.algorithms.spot_prediction import reflection_frames
    from dials.algorithms.shoebox import BBoxCalculator
    from math import pi

    # Create the spot predictor
    s0 = self.beam.get_s0()
    m2 = self.gonio.get_rotation_axis()
    UB = self.UB
    dphi = self.scan.get_oscillation_range(deg=False)
    predict_rays = RayPredictor(s0, m2, dphi)

    # Predict the reflections
    self._new_reflections = reflection_frames(self.scan, ray_intersection(
        self.detector, predict_rays(self.miller_indices, UB)))

    # Set the divergence and mosaicity
    n_sigma = 5.0
    delta_divergence = n_sigma * self.sigma_divergence * pi / 180.0
    delta_mosaicity = n_sigma * self.sigma_mosaicity * pi / 180.0

    # Create the bounding box calculator
    calculate_bbox = BBoxCalculator(self.beam, self.detector, self.gonio,
        self.scan, delta_divergence, delta_mosaicity)

    # Calculate the frame numbers of all the reflections
    calculate_bbox(self._new_reflections)

def read_reflection_file(reflection_file):
  '''Read reflections from pickle file.'''
  from dials.model.data import ReflectionList
  import cPickle as pickle
  return pickle.load(open(reflection_file, 'rb'))

def run(reflection_file):
  '''Do the refinement.'''
  from scitbx import matrix

  # Load the reflections from the pickle file
  reflections = read_reflection_file(reflection_file)



if __name__ == '__main__':
  import argparse

  parser = argparse.ArgumentParser()
  parser.add_argument("reflection_file", help="/path/to/data.pkl")
  parser.add_argument("xparm_file", help="/path/to/XPARM.XDS")
  parser.add_argument("integrate_file", help="/path/to/INTEGRATE.HKL")
  parser.add_argument("-v", "--verbosity", action="count", default=0,
                      help="set verbosity level; -vv gives verbosity level 2")
  parser.add_argument("--scan_varying",
          help="experimental option to follow refinement to " + \
               "convergence with a round of refinement using a " + \
               "scan-varying parameterisation (i.e. time dependence)",
               action="store_true")
  args = parser.parse_args()

  # reconstitute the reflections
  reflections = read_reflection_file(args.reflection_file)

  # Run the refinement
  runner = RefinementRunner(args.xparm_file, args.integrate_file,
                            reflections, args.verbosity,
                            args.scan_varying)
  runner()

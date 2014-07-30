#
#  Copyright (C) (2013) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman.
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#

"""This is deprecated code that is now used only in a regression test comparing
results to new versions of the calculations.

* ScanVaryingReflectionPredictor performs prediction for a particular image with
  different setting matrices at beginning and end of the image
* ScanVaryingReflectionListGenerator runs ScanVaryingReflectionPredictor over
  all images in a Scan

"""

from __future__ import division

from math import pi, sqrt, acos, atan2, fabs
from scitbx import matrix
import scitbx.math
from cctbx.array_family import flex

from dials.model.data import Reflection, ReflectionList

from dials.algorithms.spot_prediction.reeke import solve_quad
from dials.algorithms.spot_prediction import ReekeIndexGenerator

class ScanVaryingReflectionPredictor(object):
  """
  Reflection prediction for a relp within a small segment of a scan, which
  we assume to be a single image.

  The path of the relp through reciprocal space between the start and end of
  the image is approximated by a general linear transformation (not just a
  rotation).

  Temporarily, we initialise with a list of N+1 UB matrices, where N is the
  total number of images. In future we will simple pass a crystal model,
  which will store its own per-image UB matrix.

  Currently it is assumed that only the crystal model varies with image
  number, whilst the other models remain static.
  """

  def __init__(self, crystal, beam, gonio, scan, dmin):

    if crystal.num_scan_points == 0:
      raise TypeError("The provided crystal has no scan point samples of the setting matrix")

    self._crystal = crystal
    self._beam = beam
    self._gonio = gonio
    self._scan = scan

    # get offset to convert image array number to array index in crystal
    self._first_image = scan.get_array_range()[0]

    # resolution limit
    self._dmin = dmin
    self._dstarmax = 1. / dmin
    self._dstarmax_sq = self._dstarmax**2

    # reciprocal space beam vector
    self._s0 = matrix.col(self._beam.get_s0())
    self._s0mag_sq = self._s0.length_sq()
    self._s0mag = sqrt(self._s0mag_sq)

    # rotation axis
    self._axis = matrix.col(self._gonio.get_rotation_axis())

    # attributes to be set during preparation for a particular
    # image number
    self._step = None
    self._image_number = None

    # attributes to be set during prediction and called by user of
    # the class
    self.is_outside_resolution_limit = False

  def prepare(self, image_number, step = 1):
    """
    Cache transformations that position relps at the beginning and end of
    the step.
    """

    self._image_number = image_number
    self._step = step

    phi_beg = self._scan.get_angle_from_array_index(image_number,
                                                    deg = False)
    phi_end = self._scan.get_angle_from_array_index(image_number + step,
                                                    deg = False)
    r_beg = matrix.sqr(scitbx.math.r3_rotation_axis_and_angle_as_matrix(
        axis = self._axis, angle = phi_beg, deg = False))
    r_end = matrix.sqr(scitbx.math.r3_rotation_axis_and_angle_as_matrix(
        axis = self._axis, angle = phi_end, deg = False))

    self._A1 = r_beg * self._crystal.get_A_at_scan_point(image_number - \
                                                         self._first_image)

    self._A2 = r_end * self._crystal.get_A_at_scan_point(image_number - \
                                                      self._first_image + step)

    return

  def get_A1(self):
    """Get the setting matrix for the beginning of the step"""

    return self._A1

  def get_A2(self):
    """Get the setting matrix for the end of the step"""

    return self._A2

  def predict(self, hkl):
    """
    Predict for hkl during the passage from A1*h to A2*h.

    If a prediction is found, return the predicted Reflection. Otherwise
    return None.
    """

    self._r1 = self._A1 * matrix.col(hkl)
    self._r2 = self._A2 * matrix.col(hkl)

    dr = self._r2 - self._r1
    s0pr1 = self._s0 + self._r1

    # distances from Ewald sphere along radii
    r1_from_ES = (s0pr1).length() - self._s0mag
    r2_from_ES = (self._s0 + self._r2).length() - self._s0mag

    starts_outside_ES = r1_from_ES >= 0.
    ends_outside_ES = r2_from_ES >= 0.

    self.is_outside_res_limit = self._r1.length_sq() > self._dstarmax_sq

    # stop prediction for a relp that doesn't cross the ES (or crosses 2x)
    if starts_outside_ES == ends_outside_ES: return None

    # stop prediction if the relp is outside the resolution limit
    if self.is_outside_res_limit: return None

    # solve equation |s0 + r1 + alpha * dr| = |s0| for alpha. This is
    # equivalent to solving the quadratic equation
    #
    # alpha^2*dr.dr + 2*alpha(s0 + r1).dr + 2*s0.r1 + r1.r1 = 0

    roots = solve_quad(dr.length_sq(),
                       2*s0pr1.dot(dr),
                       self._r1.length_sq() + 2*self._s0.dot(self._r1))

    # choose a root that is in [0,1]
    alpha = filter(lambda x: x is not None and (0 <= x <= 1), roots)[0]

    # calculate the scattering vector s1
    s1 = self._r1 + alpha * dr + self._s0

    # calculate approximate frame and rotation angle
    frame = self._image_number + self._step * alpha
    angle = self._scan.get_angle_from_array_index(frame, deg = False)

    # create the Reflection and set properties
    r = Reflection(hkl)
    r.beam_vector = s1
    r.rotation_angle = angle
    r.frame_number = frame
    r.entering = starts_outside_ES

    return r


class ScanVaryingReflectionPredictorDebug(ScanVaryingReflectionPredictor):
  """Debugging version of ScanVaryingReflectionPredictor that stores all
  the sites it tests and provides the ability to write them out by abusing
  PDB format."""

  # reciprocal lattice coords of predictions
  trial_sites = []

  def predict(self, hkl):

    r = super(ScanVaryingReflectionPredictorDebug, self).predict(hkl)

    # store reciprocal lattice coordinates of the trial
    self.trial_sites.append(self._r1)

  # write out search points to PDB to view by coot
  def debug_write_reciprocal_lattice_points_as_pdb(self,
                              file_name='reciprocal_lattice.pdb'):
    from cctbx import crystal, xray
    cs = crystal.symmetry(unit_cell=(1000,1000,1000,90,90,90),
                          space_group="P1")
    xs = xray.structure(crystal_symmetry=cs)
    for site in self.trial_sites:
      xs.add_scatterer(xray.scatterer("C", site=site))
    xs.sites_mod_short()
    with open(file_name, 'wb') as f:
      print >> f, xs.as_pdb_file()

class ScanVaryingReflectionListGenerator(object):
  """
  Generate and predict all reflections for a scan with a varying crystal
  model.

  Starting from the (0,0,0) reflection, known to be on the Ewald sphere, for a
  particular image t, generate indices using the Reeke algorithm, then predict
  using ScanVaryingReflectionPredictor to test each candidate.

  Temporarily pass a UBlist (see docstring for ScanVaryingReflectionPredictor)
  until we can store per-image UB matrices in a crystal model
  """

  def __init__(self, crystal, beam,
                  gonio, scan, dmin):

    self._scan = scan
    self._s0 = matrix.col(beam.get_s0())
    self._axis = matrix.col(gonio.get_rotation_axis())
    self._dmin = dmin
    self._predictor = ScanVaryingReflectionPredictor(
                        crystal, beam, gonio, scan, dmin)
    self._reflections = []

  def __call__(self):

    # reset generated reflections list
    self._reflections = []

    # loop over images
    self.step_over_images()

    return self._reflections

  def step_over_images(self):
    """Loop over images, doing the search on each and extending the
    predictions list"""

    from libtbx import easy_mp
    from dials.util import mp


    n_images = self._scan.get_num_images()

    # Change the number of processors if necessary
    nproc = mp.nproc
    if nproc > n_images:
      nproc = n_images

    iterable = self._make_blocks(n_images, nproc)

    ref_list_of_list = easy_mp.parallel_map(
        func=self._search_on_image_range,
        iterable=iterable,
        processes=nproc,
        method=mp.method,
        preserve_order=True)

    self._reflections = [e for l in ref_list_of_list for e in l]
    return

  def _make_blocks(self, n_images, num_blocks):

    blocksizes = [n_images // num_blocks] * num_blocks

    # increase the final blocksize by the remainder
    blocksizes[-1] += n_images % num_blocks

    blockranges = []
    ar_range = self._scan.get_array_range()
    start = ar_range[0]
    for block in blocksizes:
      blockranges.append((start, start + block - 1))
      start += block

    return blockranges

  def _search_on_image_range(self, indices):

    reflections = [e for t in range(indices[0], indices[1]+1) \
                   for e in self._search_on_image(t)]
    return reflections

  def _search_on_image(self, t):

    self._predictor.prepare(t)

    A1 = self._predictor.get_A1()
    A2 = self._predictor.get_A2()

    index_generator = ReekeIndexGenerator(A1, A2, self._axis, self._s0,
                                  self._dmin, margin = 1)

    indices = index_generator.to_array()

    reflections = []
    for hkl in indices:
      r = self._predictor.predict(hkl)
      if r: reflections.append(r)
    return reflections


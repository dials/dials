#
#  Copyright (C) (2013) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman.
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#

"""Reflection prediction for refinement.

Various deprecated classes that I still want to keep a visible record of

"""

from __future__ import division

from math import pi, sqrt, acos, atan2, fabs
from scitbx import matrix
import scitbx.math
from cctbx.array_family import flex
from dials.algorithms.spot_prediction import RayPredictor
from rstbx.diffraction import reflection_prediction
from rstbx.diffraction import rotation_angles

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
    #from dials.util import mp
    n_images = self._scan.get_num_images()

    # Change the number of processors if necessary
    nproc = 1
    if nproc > n_images:
      nproc = n_images

    iterable = self._make_blocks(n_images, nproc)

    ref_list_of_list = easy_mp.parallel_map(
        func=self._search_on_image_range,
        iterable=iterable,
        processes=nproc,
        method="multiprocessing",
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
    from cctbx.sgtbx import space_group_info

    space_group_type = space_group_info("P 1").group().type()
    self._predictor.prepare(t)

    A1 = self._predictor.get_A1()
    A2 = self._predictor.get_A2()

    index_generator = ReekeIndexGenerator(A1, A2, space_group_type,
                                          self._axis, self._s0,
                                  self._dmin, margin = 1)

    indices = index_generator.to_array()

    reflections = []
    for hkl in indices:
      r = self._predictor.predict(hkl)
      if r: reflections.append(r)
    return reflections


class AnglePredictor_rstbx(object):
  """
  Predict the reflecting angles for a relp based on the current states
  of models in the experimental geometry model. This version is a wrapper
  for rstbx's C++ rotation_angles so is faster than the pure Python class
  AnglePredictor_py.

  This is essentially deprecated by DIALS reflection prediction code.
  """

  def __init__(self, crystal, beam, gonio, dmin):
    """Construct by linking to instances of experimental model classes"""

    self._crystal = crystal
    self._beam = beam
    self._gonio = gonio

    self._dmin = dmin

  # To convert from the Rossmann frame we use two of Graeme's functions
  @staticmethod
  def orthogonal_component(reference, changing):
    """Return unit vector corresponding to component of changing orthogonal
    to reference."""

    r = reference.normalize()
    c = changing.normalize()

    return (c - c.dot(r) * r).normalize()

  @staticmethod
  def align_reference_frame(primary_axis, primary_target,
                            secondary_axis, secondary_target):
    """Compute a rotation matrix R: R x primary_axis = primary_target and
    R x secondary_axis places the secondary_axis in the plane perpendicular
    to the primary_target, as close as possible to the secondary_target.
    Require: primary_target orthogonal to secondary_target, primary axis
    not colinear with secondary axis."""

    from scitbx import matrix

    if type(primary_axis) == type(()) or type(primary_axis) == type([]):
      primary_axis = matrix.col(primary_axis).normalize()
    else:
      primary_axis = primary_axis.normalize()

    if type(primary_target) == type(()) or type(primary_target) == type([]):
      primary_target = matrix.col(primary_target).normalize()
    else:
      primary_target = primary_target.normalize()

    if type(secondary_axis) == type(()) or type(secondary_axis) == type([]):
      secondary_axis = matrix.col(secondary_axis).normalize()
    else:
      secondary_axis = secondary_axis.normalize()

    if type(secondary_target) == type(()) or \
           type(secondary_target) == type([]):
      secondary_target = matrix.col(secondary_target).normalize()
    else:
      secondary_target = secondary_target.normalize()

    # check properties of input axes

    assert(fabs(primary_axis.angle(secondary_axis) % pi) > 0.001)
    assert(primary_target.dot(secondary_target) < 0.001)

    if primary_target.angle(primary_axis) % pi:
      axis_p = primary_target.cross(primary_axis)
      angle_p = - primary_target.angle(primary_axis)
      Rprimary = axis_p.axis_and_angle_as_r3_rotation_matrix(angle_p)
    elif primary_target.angle(primary_axis) < 0:
      axis_p = primary_axis.ortho().normalize()
      angle_p = pi
      Rprimary = axis_p.axis_and_angle_as_r3_rotation_matrix(angle_p)
    else:
      Rprimary = matrix.identity(3)

    Rpsa = Rprimary * secondary_axis
    axis_r = secondary_target.cross(Rpsa)
    axis_s = primary_target
    oc = AnglePredictor_rstbx.orthogonal_component
    angle_s = oc(axis_s, secondary_target).angle(oc(axis_s, Rpsa))
    if axis_r.angle(primary_target) <= 0.5 * pi:
      angle_s *= -1.

    Rsecondary = axis_s.axis_and_angle_as_r3_rotation_matrix(angle_s)

    return Rsecondary * Rprimary

  def observed_indices_and_angles_from_angle_range(self, phi_start_rad,
          phi_end_rad, indices):
    """For a list of indices, return those indices that can be rotated into
    diffracting position, along with the corresponding angles"""

    # calculate conversion matrix to rossmann frame.
    R_to_rossmann = self.align_reference_frame(
                self._beam.get_unit_s0(), (0.0, 0.0, 1.0),
                self._gonio.get_rotation_axis(), (0.0, 1.0, 0.0))

    # Create a rotation_angles object for the current geometry
    ra = rotation_angles(self._dmin,
             R_to_rossmann * self._crystal.get_U() * self._crystal.get_B(),
             self._beam.get_wavelength(),
             R_to_rossmann * matrix.col(self._gonio.get_rotation_axis()))

    (obs_indices,
     obs_angles) = ra.observed_indices_and_angles_from_angle_range(
        phi_start_rad = 0.0, phi_end_rad = pi, indices = indices)

    # convert to integer miller indices
    obs_indices_int = flex.miller_index()
    for hkl in obs_indices:
      hkl = map(lambda x: int(round(x)), hkl)
      obs_indices_int.append(tuple(hkl))
    return obs_indices_int, obs_angles

  def predict(self, hkl):
    """Solve the prediction formula for the reflecting angle phi"""

    # calculate conversion matrix to rossmann frame.
    R_to_rossmann = self.align_reference_frame(
                self._beam.get_unit_s0(), (0.0, 0.0, 1.0),
                self._gonio.get_rotation_axis(), (0.0, 1.0, 0.0))

    # Create a rotation_angles object for the current geometry
    ra = rotation_angles(self._dmin,
             R_to_rossmann * self._crystal.get_U() * self._crystal.get_B(),
             self._beam.get_wavelength(),
             R_to_rossmann * matrix.col(self._gonio.get_rotation_axis()))

    if ra(hkl):

      return ra.get_intersection_angles()

    else: return None

class AnglePredictor_py(object):
  """Predict the reflecting angles for a relp based on the current states
  of models in the experimental geometry model."""

  def __init__(self, crystal, beam, gonio, dmin):
    """Construct by linking to instances of experimental model classes"""

    self._crystal = crystal
    self._beam = beam
    self._gonio = gonio

    self._dstarmax = 1. / dmin
    self._dstarmax_sq = self._dstarmax**2


  def _prepare(self):
    """Cache required quantities that are not dependent on hkl"""

    # obtain current crystal setting matrix
    U = self._crystal.get_U()
    B = self._crystal.get_B()
    self._UB = U * B

    # obtain current reciprocal space beam vector
    self._s0 = matrix.col(self._beam.get_s0())
    self._s0mag_sq = self._s0.length_sq()
    self._s0mag = sqrt(self._s0mag_sq)

    # obtain rotation axis
    self._axis = matrix.col(self._gonio.get_rotation_axis())

    # calculate their dot product
    self._axis_dot_s0 = self._axis.dot(self._s0)

  def predict(self, hkl):
    """Solve the prediction formula for the reflecting angle phi"""

    self._prepare()

    return self._predict_core(hkl)


  def _predict_core(self, hkl):

    h = matrix.col(hkl)

    # r0 is the reciprocal lattice vector at phi = 0
    r0 = self._UB * h
    r0mag_sq = r0.length_sq()

    if r0mag_sq <= 1.e-6: return None # Relp is the recip. latt. origin
    if r0mag_sq > self._dstarmax_sq: return None # Outside resolution limit

    r0_dot_axis = r0.dot(self._axis)

    # Test if the projection of relp lies outside the Ewald sphere
    r0_dot_axis_sq = r0_dot_axis**2
    Cn_rad_sq = self._s0mag_sq - r0_dot_axis_sq
    if Cn_rad_sq <= 0.0: return None

    # Radius of circle of intersection of plane of rotation of the Relp
    # with the Ewald sphere (radius of 'Cn' cf rstbx::rotation_angles)
    Cn_rad  = sqrt(Cn_rad_sq)

    # If the greatest possible distance of the relp along the beam is larger
    # than the radius of this circle offset by the radius of the Ewald
    # sphere, then it can never intersect
    r0_beam_max = sqrt(r0mag_sq - r0_dot_axis_sq)
    if r0_beam_max > Cn_rad + self._s0mag: return None

    # Likewise, if the greatest projection of the relp along the beam is
    # smaller than the difference between the radius of the Ewald sphere and
    # the circle Cn, then it cannot intersect (too close to spindle)
    if r0_beam_max < self._s0mag - Cn_rad: return None

    r0_dot_s0 = r0.dot(self._s0)
    alpha = r0_dot_s0 - r0_dot_axis * self._axis_dot_s0
    beta = self._s0.dot(self._axis.cross(r0))
    gamma = -0.5 * r0.dot(r0) - r0_dot_axis * self._axis_dot_s0

    # The prediction formula is of form
    # alpha * cos(phi) + beta * sin(phi) = gamma. We solve this following
    # the method in DSTAR from MOSFLM
    tst = alpha**2 + beta**2
    if tst <= 0.0: return None # Relp on rotation axis
    x = gamma / sqrt(tst)
    x = min(x, 1.0)
    x = max(x, -1.0)

    t1 = acos(x)
    t2 = atan2(beta, alpha)

    # Two solutions (entrance and exit of Ewald sphere). Ensure return
    # values are positive angles
    phi1 = t2 + t1
    if phi1 <= .0: phi1 = 2 * pi + phi1
    phi2 = t2 - t1
    if phi2 <= .0: phi2 = 2 * pi + phi2
    phi = tuple(sorted([phi1, phi2]))
    return phi

  def observed_indices_and_angles_from_angle_range(self, phi_start_rad,
          phi_end_rad, indices):
    """For a list of indices, return those indices that can be rotated into
    diffracting position, along with the corresponding angles"""

    self._prepare()

    obs_ind = flex.miller_index()
    obs_ang = flex.double()

    for hkl in indices:
      angs = self._predict_core(hkl)
      if angs is None: continue

      for ang in angs:
        if ang >= phi_start_rad and ang <= phi_end_rad:
          obs_ind.append(hkl)
          obs_ang.append(ang)

    return (obs_ind, obs_ang)

class ImpactPredictor(object):
  """Predict observation position for supplied reflections and angles.

  This class is just a wrapper for RSTBX's reflection_prediction class (which
  is superseded by DIALS' ray_intersection function). A wrapper is necessary
  because reflection_prediction does not use the experimental models. This
  class keeps track of those models and instantiates a reflection_prediction
  object when required with the correct geometry.

  It is called ImpactPredictor, because ReflectionPredictor does not imply
  that the hkl is actually observed, whereas ImpactPredictor does

  """

  def __init__(self, detector, goniometer, beam, crystal):
    self._detector = detector
    self._gonio = goniometer
    self._beam = beam
    self._crystal = crystal

  def predict(self, indices, angles):

    # extract required information from the models
    sensor = self._detector.sensors()[0] # assume only one sensor for now
    axis = matrix.col(self._gonio.get_rotation_axis())
    s0 = matrix.col(self._beam.get_s0())
    UB = self._crystal.get_U() * self._crystal.get_B()

    # instantiate predictor
    rp = reflection_prediction(axis, s0, UB, sensor)

    # perform impact prediction
    # rp.predict does not like indices as miller_index. A conversion to
    # matrix.col fixes that
    temp = []
    for hkl in indices:
      hkl = matrix.col(hkl)
      temp.append(hkl)
    Hc, Xc, Yc, Phic, Sc = rp.predict(temp, angles)

    # Hc is now an an array of floating point vec3. How annoying! We want
    # integer Miller indices. Force that to be so.
    # Once we have new (fast) reflection prediction code that does not
    # return floating point indices then this whole class can be replaced
    temp = flex.miller_index()
    for hkl in Hc:
      hkl = map(lambda x: int(round(x)), hkl)
      temp.append(tuple(hkl))

    # Rescale normalised scattering vectors to the correct length
    wavelength = self._beam.get_wavelength()
    Sc = map(lambda s: matrix.col(s) / wavelength, Sc)

    return (temp, Xc, Yc, Phic, Sc)

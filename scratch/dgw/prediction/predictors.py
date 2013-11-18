#!/usr/bin/env cctbx.python

#
#  Copyright (C) (2013) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman.
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#

from __future__ import division

from math import pi, sin, cos, sqrt, acos, atan2, fabs
from scitbx import matrix
from cctbx.array_family import flex
from dials.algorithms.spot_prediction import RayPredictor
from rstbx.diffraction import reflection_prediction
from rstbx.diffraction import rotation_angles

class ReflectionPredictor(object):
  """Predict for a relp based on the current states of models in the
  experimental geometry model. This is a wrapper for DIALS' C++
  RayPredictor class, which does the real work. This class keeps track
  of the experimental geometry, and instantiates a RayPredictor when
  required.
  """

  def __init__(self, crystal, beam, gonio, sweep_range = (0, 2.*pi)):
    """Construct by linking to instances of experimental model classes"""

    self._crystal = crystal
    self._beam = beam
    self._gonio = gonio
    self._sweep_range = sweep_range
    self.update()

  def update(self):
    """Build a RayPredictor object for the current geometry"""

    self._ray_predictor = RayPredictor(self._beam.get_s0(),
                    self._gonio.get_rotation_axis(),
                    self._crystal.get_U() * self._crystal.get_B(),
                    self._sweep_range)

  def predict(self, hkl):
    """Solve the prediction formula for the reflecting angle phi"""

    return self._ray_predictor(hkl)

class AnglePredictor_rstbx(object):
  """Predict the reflecting angles for a relp based on the current states
  of models in the experimental geometry model. This version is a wrapper
  for rstbx's C++ rotation_angles so is faster than the pure Python class
  AnglePredictor_py"""

  def __init__(self, crystal, beam, gonio, dmin):
    """Construct by linking to instances of experimental model classes"""

    self._crystal = crystal
    self._beam = beam
    self._gonio = gonio

    self._dmin = dmin
    #self._dstarmax = 1. / dmin
    #self._dstarmax_sq = self._dstarmax**2

  # To convert from the Rossmann frame we use two of Graeme's
  # functions taken from use_case_xds_method/tdi.py
  def orthogonal_component(self, reference, changing):
    """Return unit vector corresponding to component of changing orthogonal to
    reference."""

    r = reference.normalize()
    c = changing.normalize()

    return (c - c.dot(r) * r).normalize()

  def align_reference_frame(self, primary_axis, primary_target,
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

    axis_r = secondary_target.cross(Rprimary * secondary_axis)
    axis_s = primary_target
    if (axis_r.angle(primary_target) > 0.5 * pi):
      angle_s = self.orthogonal_component(axis_s, secondary_target).angle(
          self.orthogonal_component(axis_s, Rprimary * secondary_axis))
    else:
      angle_s = - self.orthogonal_component(axis_s, secondary_target).angle(
          self.orthogonal_component(axis_s, Rprimary * secondary_axis))

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

    obs_indices, obs_angles = ra.observed_indices_and_angles_from_angle_range(
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

  This class is just a wrapper for RSTBX's reflection_prediction class (in
  future that class should be replaced). A wrapper is necessary because
  reflection_prediction does not use the experimental models. This class
  keeps track of those models and instantiates a reflection_prediction object
  when required with the correct geometry.

  It is called ImpactPredictor, because ReflectionPredictor does not imply
  that the hkl is actually observed, whereas ImpactPredictor does"""

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
    #Hc, Xc, Yc, Phic, Sc = rp.predict(indices, angles)
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

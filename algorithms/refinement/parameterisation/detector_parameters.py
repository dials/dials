#!/usr/bin/env python

#
#  Copyright (C) (2013) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman.
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#

from __future__ import division
from scitbx import matrix

from model_parameters import Parameter, ModelParameterisation
from dials.algorithms.refinement.refinement_helpers \
    import dR_from_axis_and_angle

class DetectorParameterisationSinglePanel(ModelParameterisation):
  """Parameterisation for a single abstract panel
  plane, with angles expressed in mrad"""

  def __init__(self, detector):

    # The state of a single Panel is its detector matrix d = (d1|d2|d0).
    # However, for the purposes of parameterisation we choose a different
    # vector than d0 to locate the Panel. That's because we want to perform
    # rotations around a point on the detector surface, and d0 points to the
    # corner of the Panel. To avoid excess correlations between 'tilt' and
    # 'twist' angles with the detector distance, we prefer to perform
    # rotations around a point located at the centre of the Panel. This is
    # usually close to point of intersection with the plane normal drawn
    # from the origin of the laboratory frame.
    #
    # Therefore we define:
    #
    # * a vector 'dorg' locating the centre of the single Panel
    # * a pair of orthogonal unit directions 'd1' and 'd2' forming a plane
    #   with its origin at the end of the vector dorg
    # * a third unit direction 'dn', orthogonal to both 'd1' & 'd2'.
    # * offsets to locate the origin d0 of the Panel frame from the
    #   tip of the dorg vector, in terms of the coordinate system
    #   formed by d1, d2 and dn.
    #
    # Held separately in attribute 'models' are:
    # * references to the detector objects contained in this model
    #
    # For this simplified class there is only a single Panel frame
    # and the vector dn is not actually required, because the plane formed
    # by d1 and d2 is coplanar with the sensor plane. Therefore the
    # offset is fully in terms of d1 and d2

    # set up the initial state of the detector parameterisation from the
    # orientation of the single Panel it contains, in terms of the vectors
    # dorg, d1 and d2.

    # get some vectors we need from the Panel
    panel = detector[0]
    so = matrix.col(panel.get_origin())
    d1 = matrix.col(panel.get_fast_axis())
    d2 = matrix.col(panel.get_slow_axis())
    dn = matrix.col(panel.get_normal())

    # We choose the dorg vector to terminate in the centre of the Panel, and
    # the offset between the end of the dorg vector and the Panel origin is
    # a coordinate matrix with elements in the basis d1, d2, dn
    panel_lim = panel.get_image_size_mm()
    offset = matrix.col((-1. * panel_lim[0] / 2.,
                         -1. * panel_lim[1] / 2.,
                          0.))

    dorg = so - offset[0] * d1 - offset[1] * d2

    # Set up the initial state. There are multiple items of interest, so
    # use a dictionary here (note for a single panel model we can do
    # without the first 3 of these, but will need them for multiple panels)
    istate = {'d1':d1,
              'd2':d2,
              'dn':dn,
              'offset':offset}

    # set up the parameters.
    # distance from lab origin to detector model plane along its normal, in
    # initial orientation
    distance = panel.get_distance()
    dist = Parameter(distance, dn, 'length (mm)', 'Dist')

    # shift in the detector model plane to locate dorg, in initial
    # orientation
    shift = dorg - dn * distance
    shift1 = Parameter(shift.dot(d1), d1, 'length (mm)', 'Shift1')
    shift2 = Parameter(shift.dot(d2), d2, 'length (mm)', 'Shift2')

    # rotations of the plane through its origin about:
    # 1) axis normal to initial orientation
    # 2) d1 axis of initial orientation
    # 3) d2 axis of initial orientation
    tau1 = Parameter(0, dn, 'angle (mrad)', 'Tau1')
    tau2 = Parameter(0, d1, 'angle (mrad)', 'Tau2')
    tau3 = Parameter(0, d2, 'angle (mrad)', 'Tau3')

    # build the parameter list in a specific,  maintained order
    p_list = [dist, shift1, shift2, tau1, tau2, tau3]

    # set up the list of model objects being parameterised (here
    # just the detector containing a single panel)
    models = [detector]

    # set up the base class
    ModelParameterisation.__init__(self, models, istate, p_list)

    # call compose to calculate all the derivatives
    self.compose()

    return

  def compose(self):

    # extract items from the initial state
    id1 = self._initial_state['d1']
    id2 = self._initial_state['d2']
    ioffset = self._initial_state['offset']

    # extract parameters from the internal list
    dist, shift1, shift2, tau1, tau2, tau3 = self._param

    # convert angles to radians
    tau1rad = tau1.value / 1000.
    tau2rad = tau2.value / 1000.
    tau3rad = tau3.value / 1000.

    # compose rotation matrices and their first order derivatives
    Tau1 = (tau1.axis).axis_and_angle_as_r3_rotation_matrix(tau1rad,
                                                            deg=False)
    dTau1_dtau1 = dR_from_axis_and_angle(tau1.axis, tau1rad, deg=False)

    Tau2 = (tau2.axis).axis_and_angle_as_r3_rotation_matrix(tau2rad,
                                                            deg=False)
    dTau2_dtau2 = dR_from_axis_and_angle(tau2.axis, tau2rad, deg=False)

    Tau3 = (tau3.axis).axis_and_angle_as_r3_rotation_matrix(tau3rad,
                                                            deg=False)
    dTau3_dtau3 = dR_from_axis_and_angle(tau3.axis, tau3rad, deg=False)

    Tau32 = Tau3 * Tau2
    Tau321 = Tau32 * Tau1

    # Compose new state
    # =================
    # First the frame positioned at a distance from the lab origin
    P0 = dist.value * dist.axis # distance along initial detector normal
    Px = P0 + id1 # point at the end of d1 in lab frame
    Py = P0 + id2 # point at the end of d2 in lab frame

    # detector shift vector
    dsv = P0 + shift1.value * shift1.axis + shift2.value * shift2.axis

    # compose dorg point
    dorg = Tau321 * dsv - Tau32 * P0 + P0

    # compose d1, d2 and dn and ensure frame remains orthonormal.
    d1 = (Tau321 * (Px - P0)).normalize()
    d2 = (Tau321 * (Py - P0)).normalize()
    dn = d1.cross(d2).normalize()
    # NB dn not actually used in this simple model; calculation left
    # here as a reminder for future extension
    d2 = dn.cross(d1)

    # compose new sensor origin
    o = dorg + ioffset[0] * d1 + ioffset[1] * d2

    # compose new panel dir1, dir2. For this model with a single
    # sensor, this is trivial as they are the same as d1 and d2
    dir1 = d1
    dir2 = d2

    # now update the panel with its new position and orientation.
    # The detector is the first model in _models, the panel is the
    # first in the detector
    (self._models[0])[0].set_frame(dir1, dir2, o)

    # calculate derivatives of the state wrt parameters
    # =================================================
    # Start with the dorg vector, where
    # dorg = Tau321 * dsv - Tau32 * P0 + P0

    # derivative wrt dist
    dP0_ddist = dist.axis
    ddsv_ddist = dP0_ddist
    ddorg_ddist = Tau321 * ddsv_ddist - Tau32 * dP0_ddist + \
                  dP0_ddist

    # derivative wrt shift1
    ddsv_dshift1 = shift1.axis
    ddorg_dshift1 = Tau321 * ddsv_dshift1

    # derivative wrt shift2
    ddsv_dshift2 = shift2.axis
    ddorg_dshift2 = Tau321 * ddsv_dshift2

    # derivative wrt tau1
    dTau321_dtau1 = Tau32 * dTau1_dtau1
    ddorg_dtau1 = dTau321_dtau1 * dsv

    # derivative wrt tau2
    dTau32_dtau2 = Tau3 * dTau2_dtau2
    dTau321_dtau2 = dTau32_dtau2 * Tau1
    ddorg_dtau2 = dTau321_dtau2 * dsv - dTau32_dtau2 * P0

    # derivative wrt tau3
    dTau32_dtau3 = dTau3_dtau3 * Tau2
    dTau321_dtau3 = dTau32_dtau3 * Tau1
    ddorg_dtau3 = dTau321_dtau3 * dsv - dTau32_dtau3 * P0

    # Now derivatives of the direction d1, where
    # d1 = (Tau321 * (Px - P0)).normalize()
    # For calc of derivatives ignore the normalize(), which should
    # be unnecessary anyway as Px - P0 is a unit vector and Tau321 a
    # pure rotation.

    # derivative wrt dist
    # dPx_ddist = dist.axis; dP0_ddist = dist.axis, so these cancel
    dd1_ddist = matrix.col((0., 0., 0.))

    # derivative wrt shift1
    dd1_dshift1 = matrix.col((0., 0., 0.))

    # derivative wrt shift2
    dd1_dshift2 = matrix.col((0., 0., 0.))

    # derivative wrt tau1
    dd1_dtau1 = dTau321_dtau1 * (Px - P0)

    # derivative wrt tau2
    dd1_dtau2 = dTau321_dtau2 * (Px - P0)

    # derivative wrt tau3
    dd1_dtau3 = dTau321_dtau3 * (Px - P0)

    # Derivatives of the direction d2, where
    # d2 = (Tau321 * (Py - P0)).normalize()

    # derivative wrt dist
    dd2_ddist = matrix.col((0., 0., 0.))

    # derivative wrt shift1
    dd2_dshift1 = matrix.col((0., 0., 0.))

    # derivative wrt shift2
    dd2_dshift2 = matrix.col((0., 0., 0.))

    # derivative wrt tau1
    dd2_dtau1 = dTau321_dtau1 * (Py - P0)

    # derivative wrt tau2
    dd2_dtau2 = dTau321_dtau2 * (Py - P0)

    # derivative wrt tau3
    dd2_dtau3 = dTau321_dtau3 * (Py - P0)

    # Derivatives of the direction dn, where
    # dn = d1.cross(d2).normalize()

    # derivative wrt dist
    ddn_ddist = matrix.col((0., 0., 0.))

    # derivative wrt shift1
    ddn_dshift1 = matrix.col((0., 0., 0.))

    # derivative wrt shift2
    ddn_dshift2 = matrix.col((0., 0., 0.))

    # derivative wrt tau1. Product rule for cross product applies
    ddn_dtau1 = dd1_dtau1.cross(d2) + d1.cross(dd2_dtau1)

    # derivative wrt tau2
    ddn_dtau2 = dd1_dtau2.cross(d2) + d1.cross(dd2_dtau2)

    # derivative wrt tau3
    ddn_dtau3 = dd1_dtau3.cross(d2) + d1.cross(dd2_dtau3)

    # calculate derivatives of the attached sensor matrix
    # ===================================================
    # sensor origin:
    # o = dorg + ioffset[0] * d1 + ioffset[1] * d2

    # derivative wrt dist
    do_ddist = ddorg_ddist + ioffset[0] * dd1_ddist + \
               ioffset[1] * dd2_ddist

    # derivative wrt shift1
    do_dshift1 = ddorg_dshift1 + ioffset[0] * dd1_dshift1 + \
                 ioffset[1] * dd2_dshift1

    # derivative wrt shift2
    do_dshift2 = ddorg_dshift2 + ioffset[0] * dd1_dshift2 + \
                 ioffset[1] * dd2_dshift2

    # derivative wrt tau1
    do_dtau1 = ddorg_dtau1 + ioffset[0] * dd1_dtau1 + \
               ioffset[1] * dd2_dtau1

    # derivative wrt tau2
    do_dtau2 = ddorg_dtau2 + ioffset[0] * dd1_dtau2 + \
               ioffset[1] * dd2_dtau2

    # derivative wrt tau3
    do_dtau3 = ddorg_dtau3 + ioffset[0] * dd1_dtau3 + \
               ioffset[1] * dd2_dtau3

    # sensor dir1:
    # dir1 = d1

    # derivative wrt dist
    ddir1_ddist = dd1_ddist

    # derivative wrt shift1
    ddir1_dshift1 = dd1_dshift1

    # derivative wrt shift2
    ddir1_dshift2 = dd1_dshift2

    # derivative wrt tau1
    ddir1_dtau1 = dd1_dtau1

    # derivative wrt tau2
    ddir1_dtau2 = dd1_dtau2

    # derivative wrt tau3
    ddir1_dtau3 = dd1_dtau3

    # sensor dir2:
    # dir2 = d2

    # derivative wrt dist
    ddir2_ddist = dd2_ddist

    # derivative wrt shift1
    ddir2_dshift1 = dd2_dshift1

    # derivative wrt shift2
    ddir2_dshift2 = dd2_dshift2

    # derivative wrt tau1
    ddir2_dtau1 = dd2_dtau1

    # derivative wrt tau2
    ddir2_dtau2 = dd2_dtau2

    # derivative wrt tau3
    ddir2_dtau3 = dd2_dtau3

    # combine these vectors together into derivatives of the sensor
    # matrix d and store them, converting angles back to mrad

    # derivative wrt dist
    self._dstate_dp[0] = matrix.sqr(ddir1_ddist.elems +
                          ddir2_ddist.elems +
                          do_ddist.elems).transpose()

    # derivative wrt shift1
    self._dstate_dp[1] = matrix.sqr(ddir1_dshift1.elems +
                          ddir2_dshift1.elems +
                          do_dshift1.elems).transpose()

    # derivative wrt shift2
    self._dstate_dp[2] = matrix.sqr(ddir1_dshift2.elems +
                          ddir2_dshift2.elems +
                          do_dshift2.elems).transpose()

    # derivative wrt tau1
    self._dstate_dp[3] = matrix.sqr(ddir1_dtau1.elems +
                          ddir2_dtau1.elems +
                          do_dtau1.elems).transpose() / 1000.

    # derivative wrt tau2
    self._dstate_dp[4] = matrix.sqr(ddir1_dtau2.elems +
                          ddir2_dtau2.elems +
                          do_dtau2.elems).transpose() / 1000.

    # derivative wrt tau3
    self._dstate_dp[5] = matrix.sqr(ddir1_dtau3.elems +
                          ddir2_dtau3.elems +
                          do_dtau3.elems).transpose() / 1000.

    return

  def get_state(self):

    # only a single panel exists, so no multi_state_elt argument is allowed
    panel = (self._models[0])[0]
    return matrix.sqr(panel.get_d_matrix())

class DetectorParameterisationMultiPanel(ModelParameterisation):
  """Experimental implementation of parameterisation for a multiple
  panel detector, treated as a single rigid block."""

  def __init__(self, detector, beam):

    # The state of each Panel in the detector model is its matrix
    # d = (d1|d2|d0). We need to define a new coordinate system rigidly
    # attached to the detector model in which to express the
    # parameterisation and compose each of the Panel states.
    #
    # We define:
    #
    # * a vector 'dorg' locating a point in laboratory space that moves with
    #   the rigid body of the detector and thus is fixed wrt each of the
    #   Panels.
    # * A pair of orthogonal unit directions 'd1' and 'd2' forming a plane
    #   with its origin at the end of the vector dorg.
    # * a third unit direction 'dn', orthogonal to both 'd1' & 'd2'.
    # * offsets to locate the origin of each panel frame from the
    #   tip of the dorg vector, in terms of the coordinate system
    #   formed by d1, d2 and dn.
    #
    # Held separately in attribute 'models' are:
    # * references to detector objects contained in this model

    # set up the initial state of the detector model from the
    # orientation of whichever Panel has its centre most closely
    # located to the direct beam intersection. Call this 'mid_panel'

    beam_centres = [matrix.col(p.get_beam_centre(beam.get_unit_s0())) \
                    for p in detector]
    panel_centres = [0.5 * matrix.col(p.get_image_size_mm())
                     for p in detector]
    beam_to_centres = [(a - b).length() for a, b in \
                      zip(beam_centres, panel_centres)]
    mid_panel_id = beam_to_centres.index(min(beam_to_centres))
    mid_panel = detector[mid_panel_id]

    # get some vectors we need from the mid_panel
    so = matrix.col(mid_panel.get_origin())
    d1 = matrix.col(mid_panel.get_fast_axis())
    d2 = matrix.col(mid_panel.get_slow_axis())
    dn = matrix.col(mid_panel.get_normal())

    # we choose the dorg vector to terminate in the centre of the mid_panel,
    # and the offset between the end of the dorg vector and each Panel
    # origin is a coordinate matrix with elements in the basis d1, d2, dn.
    # We need also each Panel's plane directions dir1 and dir2 in terms of
    # d1, d2 and dn.
    mid_panel_centre = panel_centres[mid_panel_id]
    dorg = so + mid_panel_centre[0] * d1 + mid_panel_centre[1] * d2

    offsets, dir1s, dir2s = [], [], []
    for p in detector:
      offset = matrix.col(p.get_origin()) - dorg
      offsets.append(matrix.col((offset.dot(d1),
                                 offset.dot(d2),
                                 offset.dot(dn))))
      dir1 = matrix.col(p.get_fast_axis())
      dir1_new_basis = matrix.col((dir1.dot(d1),
                                   dir1.dot(d2),
                                   dir1.dot(dn)))
      dir1s.append(dir1_new_basis)
      dir2 = matrix.col(p.get_slow_axis())
      dir2_new_basis = matrix.col((dir2.dot(d1),
                                   dir2.dot(d2),
                                   dir2.dot(dn)))
      dir2s.append(dir2_new_basis)

    # The offsets and directions in the d1, d2, dn basis are fixed
    # quantities, not dependent on parameter values.
    self._offsets = offsets
    self._dir1s = dir1s
    self._dir2s = dir2s

    # Set up the initial state. This is the basis d1, d2, dn.
    istate = {'d1':d1, 'd2':d2, 'dn':dn}

    # set up the parameters.
    # distance from lab origin to mid_panel plane along its normal,
    # in initial orientation
    distance = mid_panel.get_distance()
    dist = Parameter(distance, dn, 'length (mm)', 'Dist')

    # shift in the detector model plane to locate dorg, in initial
    # orientation
    shift = dorg - dn * distance
    shift1 = Parameter(shift.dot(d1), d1, 'length (mm)', 'Shift1')
    shift2 = Parameter(shift.dot(d2), d2, 'length (mm)', 'Shift2')

    # rotations of the plane through its origin about:
    # 1) axis normal to initial orientation
    # 2) d1 axis of initial orientation
    # 3) d2 axis of initial orientation
    tau1 = Parameter(0, dn, 'angle (mrad)', 'Tau1')
    tau2 = Parameter(0, d1, 'angle (mrad)', 'Tau2')
    tau3 = Parameter(0, d2, 'angle (mrad)', 'Tau3')

    # build the parameter list in a specific,  maintained order
    p_list = [dist, shift1, shift2, tau1, tau2, tau3]

    # set up the list of model objects being parameterised (here
    # just a single detector containing multiple panels)
    models = [detector]

    # set up the base class
    ModelParameterisation.__init__(self, models, istate, p_list,
                                   is_multi_state=True)

    # call compose to calculate all the derivatives
    self.compose()

    return

  def compose(self):

    # extract items from the initial state
    id1 = self._initial_state['d1']
    id2 = self._initial_state['d2']
    idn = self._initial_state['dn']

    # extract parameters from the internal list
    dist, shift1, shift2, tau1, tau2, tau3 = self._param

    # Extract the detector model, which is the first entry in _models.
    detector = self._models[0]

    # convert angles to radians
    tau1rad = tau1.value / 1000.
    tau2rad = tau2.value / 1000.
    tau3rad = tau3.value / 1000.

    # compose rotation matrices and their first order derivatives
    Tau1 = (tau1.axis).axis_and_angle_as_r3_rotation_matrix(tau1rad,
                                                            deg=False)
    dTau1_dtau1 = dR_from_axis_and_angle(tau1.axis, tau1rad, deg=False)

    Tau2 = (tau2.axis).axis_and_angle_as_r3_rotation_matrix(tau2rad,
                                                            deg=False)
    dTau2_dtau2 = dR_from_axis_and_angle(tau2.axis, tau2rad, deg=False)

    Tau3 = (tau3.axis).axis_and_angle_as_r3_rotation_matrix(tau3rad,
                                                            deg=False)
    dTau3_dtau3 = dR_from_axis_and_angle(tau3.axis, tau3rad, deg=False)

    Tau32 = Tau3 * Tau2
    Tau321 = Tau32 * Tau1

    # Compose new state
    # =================
    # First the frame positioned at a distance from the lab origin
    P0 = dist.value * dist.axis # distance along initial detector normal
    Px = P0 + id1 # point at the end of d1 in lab frame
    Py = P0 + id2 # point at the end of d2 in lab frame

    # detector shift vector
    dsv = P0 + shift1.value * shift1.axis + shift2.value * shift2.axis

    # compose dorg point
    dorg = Tau321 * dsv - Tau32 * P0 + P0

    # compose new d1, d2 and dn and ensure frame remains orthonormal.
    d1 = (Tau321 * (Px - P0)).normalize()
    d2 = (Tau321 * (Py - P0)).normalize()
    dn = d1.cross(d2).normalize()
    d2 = dn.cross(d1)

    # compose new Panel origins
    origins = [dorg + offset[0] * d1 + \
                      offset[1] * d2 + \
                      offset[2] * dn for offset in self._offsets]

    # compose new Panel directions
    dir1s = [vec[0] * d1 + \
             vec[1] * d2 + \
             vec[2] * dn for vec in self._dir1s]
    dir2s = [vec[0] * d1 + \
             vec[1] * d2 + \
             vec[2] * dn for vec in self._dir2s]

    # now update the panels with their new position and orientation.
    for p, dir1, dir2, org in zip(detector, dir1s, dir2s, origins):

      p.set_frame(dir1, dir2, org)

    # calculate derivatives of the state wrt parameters
    # =================================================
    # Start with the dorg vector, where
    # dorg = Tau321 * dsv - Tau32 * P0 + P0

    # derivative wrt dist
    dP0_ddist = dist.axis
    ddsv_ddist = dP0_ddist
    ddorg_ddist = Tau321 * ddsv_ddist - Tau32 * dP0_ddist + \
                  dP0_ddist

    # derivative wrt shift1
    ddsv_dshift1 = shift1.axis
    ddorg_dshift1 = Tau321 * ddsv_dshift1

    # derivative wrt shift2
    ddsv_dshift2 = shift2.axis
    ddorg_dshift2 = Tau321 * ddsv_dshift2

    # derivative wrt tau1
    dTau321_dtau1 = Tau32 * dTau1_dtau1
    ddorg_dtau1 = dTau321_dtau1 * dsv

    # derivative wrt tau2
    dTau32_dtau2 = Tau3 * dTau2_dtau2
    dTau321_dtau2 = dTau32_dtau2 * Tau1
    ddorg_dtau2 = dTau321_dtau2 * dsv - dTau32_dtau2 * P0

    # derivative wrt tau3
    dTau32_dtau3 = dTau3_dtau3 * Tau2
    dTau321_dtau3 = dTau32_dtau3 * Tau1
    ddorg_dtau3 = dTau321_dtau3 * dsv - dTau32_dtau3 * P0

    # Now derivatives of the direction d1, where
    # d1 = (Tau321 * (Px - P0)).normalize()
    # For calc of derivatives ignore the normalize(), which should
    # be unnecessary anyway as Px - P0 is a unit vector and Tau321 a
    # pure rotation.

    # derivative wrt dist
    # dPx_ddist = dist.axis; dP0_ddist = dist.axis, so these cancel
    dd1_ddist = matrix.col((0., 0., 0.))

    # derivative wrt shift1
    dd1_dshift1 = matrix.col((0., 0., 0.))

    # derivative wrt shift2
    dd1_dshift2 = matrix.col((0., 0., 0.))

    # derivative wrt tau1
    dd1_dtau1 = dTau321_dtau1 * (Px - P0)

    # derivative wrt tau2
    dd1_dtau2 = dTau321_dtau2 * (Px - P0)

    # derivative wrt tau3
    dd1_dtau3 = dTau321_dtau3 * (Px - P0)

    # Derivatives of the direction d2, where
    # d2 = (Tau321 * (Py - P0)).normalize()

    # derivative wrt dist
    dd2_ddist = matrix.col((0., 0., 0.))

    # derivative wrt shift1
    dd2_dshift1 = matrix.col((0., 0., 0.))

    # derivative wrt shift2
    dd2_dshift2 = matrix.col((0., 0., 0.))

    # derivative wrt tau1
    dd2_dtau1 = dTau321_dtau1 * (Py - P0)

    # derivative wrt tau2
    dd2_dtau2 = dTau321_dtau2 * (Py - P0)

    # derivative wrt tau3
    dd2_dtau3 = dTau321_dtau3 * (Py - P0)

    # Derivatives of the direction dn, where
    # dn = d1.cross(d2).normalize()

    # derivative wrt dist
    ddn_ddist = matrix.col((0., 0., 0.))

    # derivative wrt shift1
    ddn_dshift1 = matrix.col((0., 0., 0.))

    # derivative wrt shift2
    ddn_dshift2 = matrix.col((0., 0., 0.))

    # derivative wrt tau1. Product rule for cross product applies
    ddn_dtau1 = dd1_dtau1.cross(d2) + d1.cross(dd2_dtau1)

    # derivative wrt tau2
    ddn_dtau2 = dd1_dtau2.cross(d2) + d1.cross(dd2_dtau2)

    # derivative wrt tau3
    ddn_dtau3 = dd1_dtau3.cross(d2) + d1.cross(dd2_dtau3)

    # reset stored derivatives
    for i in range(len(self._dstate_dp)):
      self._dstate_dp[i] = [None] * len(detector)

    # calculate derivatives of the attached Panel matrices
    # ====================================================
    for panel_id, (offset, dir1_new_basis, dir2_new_basis) in enumerate(
                            zip(self._offsets, self._dir1s, self._dir2s)):

      # Panel origin:
      # o = dorg + offset[0] * d1 + offset[1] * d2 + offset[2] * dn

      # derivative wrt dist. NB only ddorg_ddist is not null! The other
      # elements are left here to aid understanding, but should be removed
      # when this class is ported to C++ for speed.
      do_ddist = ddorg_ddist + offset[0] * dd1_ddist + \
                               offset[1] * dd2_ddist + \
                               offset[2] * ddn_ddist

      # derivative wrt shift1. NB only ddorg_dshift1 is non-null.
      do_dshift1 = ddorg_dshift1 + offset[0] * dd1_dshift1 + \
                                   offset[1] * dd2_dshift1 + \
                                   offset[2] * ddn_dshift1

      # derivative wrt shift2. NB only ddorg_dshift2 is non-null.
      do_dshift2 = ddorg_dshift2 + offset[0] * dd1_dshift2 + \
                                   offset[1] * dd2_dshift2 + \
                                   offset[2] * ddn_dshift2

      # derivative wrt tau1
      do_dtau1 = ddorg_dtau1 + offset[0] * dd1_dtau1 + \
                               offset[1] * dd2_dtau1 + \
                               offset[2] * ddn_dtau1

      # derivative wrt tau2
      do_dtau2 = ddorg_dtau2 + offset[0] * dd1_dtau2 + \
                               offset[1] * dd2_dtau2 + \
                               offset[2] * ddn_dtau2

      # derivative wrt tau3
      do_dtau3 = ddorg_dtau3 + offset[0] * dd1_dtau3 + \
                               offset[1] * dd2_dtau3 + \
                               offset[2] * ddn_dtau3

      # Panel dir1:
      # dir1 = dir1_new_basis[0] * d1 + dir1_new_basis[1] * d2 +
      #        dir1_new_basis[2] * dn

      # derivative wrt dist. NB These are all null.
      ddir1_ddist = dir1_new_basis[0] * dd1_ddist + \
                    dir1_new_basis[1] * dd2_ddist + \
                    dir1_new_basis[2] * ddn_ddist

      # derivative wrt shift1. NB These are all null.
      ddir1_dshift1 = dir1_new_basis[0] * dd1_dshift1 + \
                      dir1_new_basis[1] * dd2_dshift1 + \
                      dir1_new_basis[2] * ddn_dshift1

      # derivative wrt shift2. NB These are all null.
      ddir1_dshift2 = dir1_new_basis[0] * dd1_dshift2 + \
                      dir1_new_basis[1] * dd2_dshift2 + \
                      dir1_new_basis[2] * ddn_dshift2

      # derivative wrt tau1
      ddir1_dtau1 = dir1_new_basis[0] * dd1_dtau1 + \
                    dir1_new_basis[1] * dd2_dtau1 + \
                    dir1_new_basis[2] * ddn_dtau1

      # derivative wrt tau2
      ddir1_dtau2 = dir1_new_basis[0] * dd1_dtau2 + \
                    dir1_new_basis[1] * dd2_dtau2 + \
                    dir1_new_basis[2] * ddn_dtau2

      # derivative wrt tau3
      ddir1_dtau3 = dir1_new_basis[0] * dd1_dtau3 + \
                    dir1_new_basis[1] * dd2_dtau3 + \
                    dir1_new_basis[2] * ddn_dtau3

      # Panel dir2:
      # dir2 = dir2_new_basis[0] * d1 + dir2_new_basis[1] * d2 +
      #        dir2_new_basis[2] * dn

      # derivative wrt dist. NB These are all null.
      ddir2_ddist = dir2_new_basis[0] * dd1_ddist + \
                    dir2_new_basis[1] * dd2_ddist + \
                    dir2_new_basis[2] * ddn_ddist

      # derivative wrt shift1. NB These are all null.
      ddir2_dshift1 = dir2_new_basis[0] * dd1_dshift1 + \
                      dir2_new_basis[1] * dd2_dshift1 + \
                      dir2_new_basis[2] * ddn_dshift1

      # derivative wrt shift2. NB These are all null.
      ddir2_dshift2 = dir2_new_basis[0] * dd1_dshift2 + \
                      dir2_new_basis[1] * dd2_dshift2 + \
                      dir2_new_basis[2] * ddn_dshift2

      # derivative wrt tau1
      ddir2_dtau1 = dir2_new_basis[0] * dd1_dtau1 + \
                    dir2_new_basis[1] * dd2_dtau1 + \
                    dir2_new_basis[2] * ddn_dtau1

      # derivative wrt tau2
      ddir2_dtau2 = dir2_new_basis[0] * dd1_dtau2 + \
                    dir2_new_basis[1] * dd2_dtau2 + \
                    dir2_new_basis[2] * ddn_dtau2

      # derivative wrt tau3
      ddir2_dtau3 = dir2_new_basis[0] * dd1_dtau3 + \
                    dir2_new_basis[1] * dd2_dtau3 + \
                    dir2_new_basis[2] * ddn_dtau3

      # combine these vectors together into derivatives of the panel
      # matrix d and store them, converting angles back to mrad

      # derivative wrt dist
      self._dstate_dp[0][panel_id] = matrix.sqr(ddir1_ddist.elems +
                  ddir2_ddist.elems + do_ddist.elems).transpose()

      # derivative wrt shift1
      self._dstate_dp[1][panel_id] = matrix.sqr(ddir1_dshift1.elems +
                  ddir2_dshift1.elems + do_dshift1.elems).transpose()

      # derivative wrt shift2
      self._dstate_dp[2][panel_id] = matrix.sqr(ddir1_dshift2.elems +
                  ddir2_dshift2.elems + do_dshift2.elems).transpose()

      # derivative wrt tau1
      self._dstate_dp[3][panel_id] = matrix.sqr(ddir1_dtau1.elems +
                  ddir2_dtau1.elems + do_dtau1.elems).transpose() / 1000.

      # derivative wrt tau2
      self._dstate_dp[4][panel_id] = matrix.sqr(ddir1_dtau2.elems +
                  ddir2_dtau2.elems + do_dtau2.elems).transpose() / 1000.

      # derivative wrt tau3
      self._dstate_dp[5][panel_id] = matrix.sqr(ddir1_dtau3.elems +
                  ddir2_dtau3.elems + do_dtau3.elems).transpose() / 1000.

    return

  def get_state(self, multi_state_elt=0):

    # There is only one detector, but the req. panel must be specified
    panel = (self._models[0])[multi_state_elt]
    return matrix.sqr(panel.get_d_matrix())

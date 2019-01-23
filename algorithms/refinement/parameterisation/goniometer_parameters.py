#!/usr/bin/env cctbx.python

#
#  Copyright (C) (2017) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman.
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#

from __future__ import absolute_import, division
from dials.algorithms.refinement.parameterisation.model_parameters import Parameter, ModelParameterisation
from scitbx import matrix
from dials.algorithms.refinement.refinement_helpers \
    import dR_from_axis_and_angle

class GoniometerMixin(object):
  """Mix-in class defining some functionality unique to goniometer
   parameterisations that can be shared by static and scan-varying versions"""

  @staticmethod
  def _build_p_list(axis, beam, parameter_type=Parameter):
    """Build the list of parameters, using the parameter_type callback to
    select between versions of the Parameter class"""

    # Set up the parameters
    if beam:
      s0u = matrix.col(beam.get_unit_s0())
      # define axis orthogonal to spindle-beam plane
      dir2 = axis.cross(s0u).normalize()
      # define axis in the spindle-beam plane
      dir1 = dir2.cross(axis).normalize()
    else:
      dir1 = axis.ortho().normalize()
      dir2 = axis.cross(dir1).normalize()

    # rotation around dir1
    gamma1 = parameter_type(.0, dir1, 'angle (mrad)', 'Gamma1')
    # rotation around dir2
    gamma2 = parameter_type(.0, dir2, 'angle (mrad)', 'Gamma2')

    # build the parameter list in a specific,  maintained order
    p_list = [gamma1, gamma2]

    return p_list

  @staticmethod
  def _compose_core(iS, gamma1, gamma2, gamma1_axis, gamma2_axis):

    # convert angles to radians
    g1rad, g2rad = gamma1 / 1000., gamma2 / 1000.

    # compose rotation matrices and their first order derivatives
    G1 = (gamma1_axis).axis_and_angle_as_r3_rotation_matrix(g1rad, deg=False)
    dG1_dg1 = dR_from_axis_and_angle(gamma1_axis, g1rad, deg=False)

    G2 = (gamma2_axis).axis_and_angle_as_r3_rotation_matrix(g2rad, deg=False)
    dG2_dg2 = dR_from_axis_and_angle(gamma2_axis, g2rad, deg=False)

    # compose new state
    G21 = G2 * G1
    S = G21 * iS

    # calculate derivatives of S wrt angles:
    #  1) derivative wrt gamma1
    dG21_dg1 = G2 * dG1_dg1
    dS_dg1 = dG21_dg1 * iS

    #  2) derivative wrt gamma2
    dG21_dg2 = dG2_dg2 * G1
    dS_dg2 = dG21_dg2 * iS

    # Convert derivatives back to mrad
    dS_dval = [dS_dg1 / 1000., dS_dg2 / 1000.]

    return S, dS_dval

class GoniometerParameterisation(ModelParameterisation, GoniometerMixin):
  """A parameterisation of a Goniometer model's setting rotation.

  The setting rotation matrix is parameterised using two orientation angles
  expressed in mrad. A beam can be provided (if present in the experiment) to
  ensure a consistent definition of the orientation angles with respect to the
  initial spindle-beam plane in the laboratory frame"""

  def __init__(self, goniometer, beam=None, experiment_ids=None):
    """Initialise the GoniometerParameterisation object

    Args:
        goniometer: A dxtbx Beam object to be parameterised.
        beam: An optional dxtbx Beam object. Defaults to None.
        experiment_ids (list): The experiment IDs affected by this
            parameterisation. Defaults to None, which is replaced by [0].
    """
    # The state of the goniometer model consists of the setting matrix [S] that
    # determines the orientation of the rotation axis in the laboratory frame
    # e_lab = [S] e_datum.
    # This parameters are orientation angles around two axes orthogonal to the
    # initial direction of e_lab. If a beam is supplied, one of these axes will
    # be within the plane containing s0 and e_lab.

    # Set up the initial state
    if experiment_ids is None:
      experiment_ids = [0]
    e_lab = matrix.col(goniometer.get_rotation_axis())
    istate = matrix.sqr(goniometer.get_setting_rotation())

    # build the parameter list
    p_list = self._build_p_list(e_lab, beam)

    # set up the base class
    ModelParameterisation.__init__(self, goniometer, istate, p_list,
                                   experiment_ids=experiment_ids)

    # call compose to calculate all the derivatives
    self.compose()

    return

  def compose(self):

    # extract setting matrix from the initial state
    iS = self._initial_state

    # extract parameters from the internal list
    gamma1, gamma2 = self._param

    # calculate new [S] and derivatives
    S, self._dstate_dp = self._compose_core(iS, gamma1.value, gamma2.value,
      gamma1_axis=gamma1.axis, gamma2_axis=gamma2.axis)

    # now update the model with its new [S]
    self._model.set_setting_rotation(S)

    return

  def get_state(self):

    # only a single setting matrix exists, so no multi_state_elt argument is
    # allowed

    return matrix.sqr(self._model.get_setting_rotation())

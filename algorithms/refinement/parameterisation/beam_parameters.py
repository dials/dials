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
from model_parameters import Parameter, ModelParameterisation
from scitbx import matrix
from dials.algorithms.refinement.refinement_helpers \
    import dR_from_axis_and_angle

class BeamParameterisation(ModelParameterisation):
  """Implementation of parameterisation for the beam with angles expressed in
  mrad and wavenumber in inverse Angstroms.

  Pass in a goniometer (if present) to ensure consistent definition of the
  beam rotation angles with respect to the spindle-beam plane."""

  def __init__(self, beam, goniometer=None, experiment_ids=None):

    # The state of the beam model consists of the s0 vector that it is
    # modelling. The initial state is the direction of this vector at the point
    # of initialisation. Future states are composed by rotations around axes
    # perpendicular to that direction and normalisation specified by the
    # wavenumber (inverse wavelength).
    #
    # The 'models' attribute refers to the beam vector contained by this
    # model.

    ### Set up the initial state
    if experiment_ids is None:
      experiment_ids = [0]
    s0 = matrix.col(beam.get_s0())
    s0dir = matrix.col(beam.get_unit_s0())
    istate = s0dir

    ### Set up the parameters
    if goniometer:
      spindle = matrix.col(goniometer.get_rotation_axis())
      s0_plane_dir2 = s0.cross(spindle).normalize()
      s0_plane_dir1 = s0_plane_dir2.cross(s0).normalize()
    else:
      s0_plane_dir1 = s0.ortho().normalize()
      s0_plane_dir2 = s0.cross(s0_plane_dir1).normalize()

    # rotation around s0_plane_dir1
    mu1 = Parameter(.0, s0_plane_dir1, 'angle (mrad)', 'Mu1')
    # rotation around s0_plane_dir2
    mu2 = Parameter(.0, s0_plane_dir2, 'angle (mrad)', 'Mu2')
    # length of s0
    nu = Parameter(s0.length(), ptype='wavenumber (Angstroem^-1)', name='nu')

    # build the parameter list in a specific,  maintained order
    p_list = [mu1, mu2, nu]

    # set up the base class
    ModelParameterisation.__init__(self, beam, istate, p_list,
                                   experiment_ids=experiment_ids)

    # call compose to calculate all the derivatives
    self.compose()

    return

  def compose(self):

    # extract direction from the initial state
    is0 = self._initial_state

    # extract parameters from the internal list
    mu1, mu2, nu = self._param

    # convert angles to radians
    mu1rad, mu2rad = mu1.value / 1000., mu2.value / 1000.

    # compose rotation matrices and their first order derivatives
    Mu1 = (mu1.axis).axis_and_angle_as_r3_rotation_matrix(mu1rad, deg=False)
    dMu1_dmu1 = dR_from_axis_and_angle(mu1.axis, mu1rad, deg=False)

    Mu2 = (mu2.axis).axis_and_angle_as_r3_rotation_matrix(mu2rad, deg=False)
    dMu2_dmu2 = dR_from_axis_and_angle(mu2.axis, mu2rad, deg=False)

    Mu21 = Mu2 * Mu1

    ### Compose new state
    s0_new_dir = (Mu21 * is0).normalize()

    # now update the model with its new s0
    self._model.set_s0(nu.value * s0_new_dir)

    ### calculate derivatives of the beam direction wrt angles
    # derivative wrt mu1
    dMu21_dmu1 = Mu2 * dMu1_dmu1
    ds0_new_dir_dmu1 = dMu21_dmu1 * is0

    # derivative wrt mu2
    dMu21_dmu2 = dMu2_dmu2 * Mu1
    ds0_new_dir_dmu2 = dMu21_dmu2 * is0

    ### calculate derivatives of the attached beam vector, converting
    ### parameters back to mrad, and store
    # derivative wrt mu1
    self._dstate_dp[0] = ds0_new_dir_dmu1 * nu.value / 1000.
    # derivative wrt mu2
    self._dstate_dp[1] = ds0_new_dir_dmu2 * nu.value / 1000.
    # derivative wrt nu
    self._dstate_dp[2] = s0_new_dir

    return

  def get_state(self):

    # only a single beam exists, so no multi_state_elt argument is allowed

    return matrix.col(self._model.get_s0())

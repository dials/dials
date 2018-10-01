#!/usr/bin/env cctbx.python

#
#  Copyright (C) (2013) STFC Rutherford Appleton Laboratory, UK.
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

class BeamMixin(object):
  """Mix-in class defining some functionality unique to beam parameterisations
  that can be shared by static and scan-varying versions"""

  @staticmethod
  def _build_p_list(s0, goniometer, parameter_type=Parameter):
    """Build the list of parameters, using the parameter_type callback to
    select between versions of the Parameter class"""

    # Set up the parameters
    if goniometer:
      spindle = matrix.col(goniometer.get_rotation_axis())
      s0_plane_dir2 = s0.cross(spindle).normalize()
      s0_plane_dir1 = s0_plane_dir2.cross(s0).normalize()
    else:
      s0_plane_dir1 = s0.ortho().normalize()
      s0_plane_dir2 = s0.cross(s0_plane_dir1).normalize()

    # rotation around s0_plane_dir1
    mu1 = parameter_type(.0, s0_plane_dir1, 'angle (mrad)', 'Mu1')
    # rotation around s0_plane_dir2
    mu2 = parameter_type(.0, s0_plane_dir2, 'angle (mrad)', 'Mu2')
    # length of s0
    nu = parameter_type(s0.length(), axis=None,
      ptype='wavenumber (Angstroem^-1)', name='nu')

    # build the parameter list in a specific,  maintained order
    p_list = [mu1, mu2, nu]

    return p_list

  @staticmethod
  def _compose_core(is0, mu1, mu2, nu, mu1_axis, mu2_axis):

    # convert angles to radians
    mu1rad, mu2rad = mu1 / 1000., mu2 / 1000.

    # compose rotation matrices and their first order derivatives
    Mu1 = (mu1_axis).axis_and_angle_as_r3_rotation_matrix(mu1rad, deg=False)
    dMu1_dmu1 = dR_from_axis_and_angle(mu1_axis, mu1rad, deg=False)

    Mu2 = (mu2_axis).axis_and_angle_as_r3_rotation_matrix(mu2rad, deg=False)
    dMu2_dmu2 = dR_from_axis_and_angle(mu2_axis, mu2rad, deg=False)

    # compose new state
    Mu21 = Mu2 * Mu1
    s0_new_dir = (Mu21 * is0).normalize()
    s0 = nu * s0_new_dir

    # calculate derivatives of the beam direction wrt angles:
    #  1) derivative wrt mu1
    dMu21_dmu1 = Mu2 * dMu1_dmu1
    ds0_new_dir_dmu1 = dMu21_dmu1 * is0

    #  2) derivative wrt mu2
    dMu21_dmu2 = dMu2_dmu2 * Mu1
    ds0_new_dir_dmu2 = dMu21_dmu2 * is0

    # calculate derivatives of the attached beam vector, converting
    # parameters back to mrad
    ds0_dval = [ds0_new_dir_dmu1 * nu / 1000.,
                ds0_new_dir_dmu2 * nu / 1000.,
                s0_new_dir]

    return s0, ds0_dval

class BeamParameterisation(ModelParameterisation, BeamMixin):
  """A parameterisation of a Beam model.

  The Beam direction and energy are parameterised using angles expressed in
  mrad and wavenumber in inverse Angstroms. A goniometer can be provided (if
  present in the experiment) to ensure a consistent definition of the beam
  rotation angles with respect to the spindle-beam plane."""

  def __init__(self, beam, goniometer=None, experiment_ids=None):
    """Initialise the BeamParameterisation object

    Args:
        beam: A dxtbx Beam object to be parameterised.
        goniometer: An optional dxtbx Goniometer object. Defaults to None.
        experiment_ids (list): The experiment IDs affected by this
            parameterisation. Defaults to None, which is replaced by [0].
    """
    # The state of the beam model consists of the s0 vector that it is
    # modelling. The initial state is the direction of this vector at the point
    # of initialisation. Future states are composed by rotations around axes
    # perpendicular to that direction and normalisation specified by the
    # wavenumber (inverse wavelength).

    # Set up the initial state
    if experiment_ids is None:
      experiment_ids = [0]
    s0 = matrix.col(beam.get_s0())
    s0dir = matrix.col(beam.get_unit_s0())
    istate = s0dir

    # build the parameter list
    p_list = self._build_p_list(s0, goniometer)

    # set up the base class
    ModelParameterisation.__init__(self, beam, istate, p_list,
                                   experiment_ids=experiment_ids)

    # call compose to calculate all the derivatives
    self.compose()

    return

  def compose(self):
    """See base class."""

    # extract direction from the initial state
    is0 = self._initial_state

    # extract parameters from the internal list
    mu1, mu2, nu = self._param

    # calculate new s0 and derivatives
    s0, self._dstate_dp = self._compose_core(is0, mu1.value, mu2.value,
      nu.value, mu1_axis=mu1.axis, mu2_axis=mu2.axis)

    # now update the model with its new s0
    self._model.set_s0(s0)

    return

  def get_state(self):
    """See base class."""

    # only a single beam exists, so no multi_state_elt argument is allowed
    return matrix.col(self._model.get_s0())

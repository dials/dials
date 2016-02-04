#!/usr/bin/env python
#
#  restraints.py
#
#  Copyright (C) 2015 Diamond Light Source and STFC Rutherford Appleton
#                     Laboratory, UK.
#
#  Author: David Waterman
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division
from scitbx.array_family import flex
from math import pi, sin, cos, sqrt
from scitbx.math import angle_derivative_wrt_vectors
from dials_refinement_helpers_ext import CalculateCellGradients
from logging import warning
DEG2RAD = pi/180.0
RAD2DEG = 180.0/pi

class DerivedParameterTie(object):
  """Calculate the restraint and gradients for a single derived parameter
  of the model"""

  def __init__(self, target, weight):

    self._target = target
    self._w = weight
    self._dRdp = None

    return

  def residual(self, parameter_value, parameter_gradients):
    '''Calculate residual R, cache gradients'''

    d = parameter_value - self._target
    self._dRdp = parameter_gradients

    return d

  def gradient(self):
    """Return dR/dp"""
    return self._dRdp

  def weight(self):
    '''Return restraint weight'''

    return self._w

class SingleUnitCellTie(object):
  """Tie the parameters of a single unit cell model parameterisation to
  target values via least-squares restraints. The restraints will be expressed
  in terms of real space unit cell constants, whilst the underlying parameters
  are encapsulated in the model parameterisation object"""

  def __init__(self, model_parameterisation, target, sigma):
    """model_parameterisation is a CrystalUnitCellParameterisation

    target is a sequence of 6 elements describing the target cell parameters

    sigma is a sequence of 6 elements giving the 'sigma' for each of the
    terms in target, from which weights for the residuals will be calculated.
    Values of zero  will remove the restraint for the cell parameter at
    that position"""

    self._xlucp = model_parameterisation
    self._target = target

    assert len(self._target) == 6
    assert len(sigma) == 6

    # calculate gradients of cell parameters wrt model parameters.
    grads = self._calculate_uc_gradients()

    # identify cell dimensions constrained to be equal
    a, b, c, aa, bb, cc = self._xlucp.get_model().get_unit_cell().parameters()
    if abs(a - b) < 1e-10:
      grad_diff = [abs(e1 - e2) for (e1, e2) in zip(grads[0], grads[1])]
      if max(grad_diff) < 1e-10:
        # a and b are equal, therefore keep only the strongest restraint
        strong, weak = sorted([sigma[0], sigma[1]])
        if strong == 0.0: strong = weak
        sigma[0] = strong
        sigma[1] = 0.0
    if abs(a - c) < 1e-10:
      grad_diff = [abs(e1 - e2) for (e1, e2) in zip(grads[0], grads[2])]
      if max(grad_diff) < 1e-10:
        # a and c are equal, therefore keep only the strongest restraint
        strong, weak = sorted([sigma[0], sigma[2]])
        if strong == 0.0: strong = weak
        sigma[0] = strong
        sigma[2] = 0.0
    if abs(b - c) < 1e-10:
      grad_diff = [abs(e1 - e2) for (e1, e2) in zip(grads[1], grads[2])]
      if max(grad_diff) < 1e-10:
        # b and c are equal, therefore keep only the strongest restraint
        strong, weak = sorted([sigma[1], sigma[2]])
        if strong == 0.0: strong = weak
        sigma[1] = strong
        sigma[2] = 0.0

    # A gradient of zero indicates that cell parameter is constrained and thus
    # to be ignored in restraints
    _sigma = []
    for sig, grad in zip(sigma, grads):
      tst = [abs(g) > 1.e-10 for g in grad]
      if any(tst):
        if sig == 0.0: sig = None
        _sigma.append(sig)
      else:
        _sigma.append(None)

    # For each non-zero sigma create a restraint between the relevant cell
    # parameter and its target value
    self._ties = []
    for t, s  in zip(self._target, _sigma):
      if s is not None:
        self._ties.append(DerivedParameterTie(t, 1./s**2))
      else:
        self._ties.append(None)

    # set up empty weights list
    self._weights = []

    return

  def _calculate_uc_gradients(self, sel=[True]*6):
    '''Calculate gradients of the unit cell parameters with respect to
    each of the parameters of the crystal unit cell model parameterisation'''

    B = self._xlucp.get_state()
    dB_dp = flex.mat3_double(self._xlucp.get_ds_dp())

    # Use C++ function for speed
    ccg = CalculateCellGradients(B, dB_dp)

    nparam = len(dB_dp)
    da = list(ccg.da_dp()) if sel[0] else [0.0] * nparam
    db = list(ccg.db_dp()) if sel[1] else [0.0] * nparam
    dc = list(ccg.dc_dp()) if sel[2] else [0.0] * nparam
    daa = list(ccg.daa_dp()) if sel[3] else [0.0] * nparam
    dbb = list(ccg.dbb_dp()) if sel[4] else [0.0] * nparam
    dcc = list(ccg.dcc_dp()) if sel[5] else [0.0] * nparam

    return (da, db, dc, daa, dbb, dcc)

  def residuals(self):
    """Calculate and return the residuals, cache gradients"""

    cell_params = self._xlucp.get_model().get_unit_cell().parameters()

    # gradients of the cell parameters wrt model parameters
    grads = self._calculate_uc_gradients(sel=[t is not None for t in self._ties])

    R = []
    for p, g, t in zip(cell_params, grads, self._ties):
      if t is None: continue
      R.append(t.residual(parameter_value=p, parameter_gradients=g))

    return R

  def gradients(self):
    """For each residual, return the gradients dR/dp. Requires residuals to be
    called first"""

    dRdp = []
    for t in self._ties:
      if t is None: continue
      dRdp.append(t.gradient())

    return dRdp

  def weights(self):
    '''Return the weights for the residuals vector'''

    # the weights do not change so cache them
    if not self._weights:
      self._weights = []
      for t in self._ties:
        if t is None: continue
        self._weights.append(t.weight())

    return self._weights

class MeanUnitCellTie(object):
  """Tie the parameters of multiple unit cell model parameterisations to
  central values via least-squares restraints. The restraints will be expressed
  in terms of real space unit cell constants, whilst the underlying parameters
  are encapsulated in the model parameterisation objects"""

  def __init__(self, model_parameterisations, sigma):
    """model_parameterisations is a list of CrystalUnitCellParameterisations

    sigma is a sequence of 6 elements giving the 'sigma' for each of the
    unit cell parameters, from which weights for the residuals will be
    calculated. Values of zero in sigma will remove the restraint for the
    cell parameter at that position"""

    self._xlucp = model_parameterisations
    self._nxls = len(model_parameterisations)

    # common factors used in gradient calculations
    self._meangradfac = 1./self._nxls
    self._gradfac = (1. - self._meangradfac)

    self._weights = []

    # initially want to calculate all gradients
    self._sel = [True] * 6

    # identify any cell dimensions constrained to be equal. If any are and a
    # restraint has been requested for that cell dimension, remove the restraint
    # for all crystals and warn in the log
    msg = ('Unit cell similarity restraints were requested for both the '
           '{0} and {1} dimensions, however for the crystal in experiment '
           '{2} these are constrained to be equal. Only the strongest '
           'of these restraints will be retained for all crystals in '
           'the restrained group.')
    for xlucp, grads in zip(self._xlucp, self.gradients()):
      a, b, c, aa, bb, cc = xlucp.get_model().get_unit_cell().parameters()
      if abs(a - b) < 1e-10:
        grad_diff = [abs(e1 - e2) for (e1, e2) in zip(grads[0], grads[1])]
        if max(grad_diff) < 1e-10:
          # a and b are equal for this crystal, therefore keep only the
          # strongest requested restraint
          if sigma[0] > 0.0 and sigma[1] > 0.0:
            warning(msg.format('a', 'b', xlucp.get_experiment_ids()[0]))
            strong, weak = sorted([sigma[0], sigma[1]])
            sigma[0] = strong
            sigma[1] = 0.0
      if abs(a - c) < 1e-10:
        grad_diff = [abs(e1 - e2) for (e1, e2) in zip(grads[0], grads[2])]
        if max(grad_diff) < 1e-10:
          # a and c are equal for this crystal, therefore keep only the
          # strongest requested restraint
          if sigma[0] > 0.0 and sigma[2] > 0.0:
            warning(msg.format('a', 'c', xlucp.get_experiment_ids()[0]))
            strong, weak = sorted([sigma[0], sigma[2]])
            sigma[0] = strong
            sigma[2] = 0.0
      if abs(b - c) < 1e-10:
        grad_diff = [abs(e1 - e2) for (e1, e2) in zip(grads[1], grads[2])]
        if max(grad_diff) < 1e-10:
          # b and c are equal for this crystal, therefore keep only the
          # strongest requested restraint
          if sigma[1] > 0.0 and sigma[2] > 0.0:
            strong, weak = sorted([sigma[1], sigma[2]])
            sigma[1] = strong
            sigma[2] = 0.0

      # A gradient of zero indicates that cell parameter is constrained and thus
      # to be ignored in restraints
      #_sigma = []
      msg = ('Unit cell similarity restraints were requested for the {0} '
             'parameter, however for the crystal in experiment {1}, {0} is '
             'constrained. This restraint will be removed for all crystals in '
             'the restrained group.')
      for i, (grad, pname) in enumerate(zip(grads, ['a', 'b', 'c', 'alpha', 'beta', 'gamma'])):
        tst = [abs(g) <= 1.e-10 for g in grad]
        if all(tst):
          # this parameter is constrained, so remove any requested restraints
          # at this position
          if sigma[i] > 0.0:
            warning(msg.format(xlucp.get_experiment_ids()[0], pname[i]))
            sigma[i] = 0.0

    # set the selection for gradient calculations to the unconstrained parameters
    self._sel = [s > 0.0 for s in sigma]

    self.nrestraints_per_cell = self._sel.count(True)

    # repeat the weights for each unit cell being restrained
    weights = flex.double([1./s**2 for s in sigma if s > 0.0])
    self._weights = flex.double()
    for xlucp in self._xlucp: self._weights.extend(weights)

    return

  def residuals(self):
    """Calculate and return the residuals"""

    cells = [xlucp.get_model().get_unit_cell().parameters() for xlucp in self._xlucp]
    a, b, c, aa, bb, cc = [flex.double(e) for e in zip(*cells)]
    resid_a = a - flex.mean(a) if self._sel[0] else None
    resid_b = b - flex.mean(b) if self._sel[1] else None
    resid_c = c - flex.mean(c) if self._sel[2] else None
    resid_aa = aa - flex.mean(aa) if self._sel[3] else None
    resid_bb = bb - flex.mean(bb) if self._sel[4] else None
    resid_cc = cc - flex.mean(cc) if self._sel[5] else None

    # collect the residuals for restrained parameters only
    resid = [e for e in [resid_a, resid_b, resid_c,
                         resid_aa, resid_bb, resid_cc] if e is not None]

    # flex array trickery to interlace the residuals
    ncol = len(resid)
    R = flex.double(flex.grid(self._nxls, ncol))
    for i, r in enumerate(resid):
      R.matrix_paste_column_in_place(r, i)

    return R.as_1d()

  def gradients(self):
    """A generator function to return the gradients dR/dp for all the restraints
    referring to a particular crystal's cell parameters.
    This only returns the gradients with respect to a single crystal
    (the one referred to by that residual). Other residuals have non-zero
    gradients however, due to the cell parameter target values being the mean
    over all cells. The gradients of the mean parameter are also set and can
    be accessed with a separate method"""

    for xlucp in self._xlucp:
      B = xlucp.get_state()
      dB_dp = flex.mat3_double(xlucp.get_ds_dp())
      # Use C++ function for speed
      ccg = CalculateCellGradients(B, dB_dp)
      dRdp = []
      self._dmeandp = []
      if self._sel[0]:
        tmp = ccg.da_dp()
        dRdp.append(list(tmp * self._gradfac))
      if self._sel[1]:
        tmp = ccg.db_dp()
        dRdp.append(list(tmp * self._gradfac))
      if self._sel[2]:
        tmp = ccg.dc_dp()
        dRdp.append(list(tmp * self._gradfac))
      if self._sel[3]:
        tmp = ccg.daa_dp()
        dRdp.append(list(tmp * self._gradfac))
      if self._sel[4]:
        tmp = ccg.dbb_dp()
        dRdp.append(list(tmp * self._gradfac))
      if self._sel[5]:
        tmp = ccg.dcc_dp()
        dRdp.append(list(tmp * self._gradfac))

      yield dRdp

  def gradients_of_the_mean(self, icell):
    """Return the gradients of the mean values of cell parameters for all the
    restraints referring to a particular unit cell."""

    xlucp = self._xlucp[icell]
    B = xlucp.get_state()
    dB_dp = flex.mat3_double(xlucp.get_ds_dp())
    # Use C++ function for speed
    ccg = CalculateCellGradients(B, dB_dp)
    dmeandp = []
    if self._sel[0]:
      tmp = ccg.da_dp()
      dmeandp.append(list(tmp * self._meangradfac))
    if self._sel[1]:
      tmp = ccg.db_dp()
      dmeandp.append(list(tmp * self._meangradfac))
    if self._sel[2]:
      tmp = ccg.dc_dp()
      dmeandp.append(list(tmp * self._meangradfac))
    if self._sel[3]:
      tmp = ccg.daa_dp()
      dRdp.append(list(tmp * self._gradfac))
      dmeandp.append(list(tmp * self._meangradfac))
    if self._sel[4]:
      tmp = ccg.dbb_dp()
      dmeandp.append(list(tmp * self._meangradfac))
    if self._sel[5]:
      tmp = ccg.dcc_dp()
      dmeandp.append(list(tmp * self._meangradfac))

    return dmeandp

  def weights(self):
    '''Return the weights for the residuals vector'''

    return self._weights

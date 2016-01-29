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
    Values of zero or None in sigma will remove the restraint for the cell
    parameter at that position"""

    self._xlucp = model_parameterisation
    self._target = target

    assert len(self._target) == 6
    assert len(sigma) == 6

    # calculate gradients of cell parameters wrt model parameters. A gradient
    # of zero means that cell parameter is constrained and thus ignored in
    # restraints
    grads = self._calculate_uc_gradients()
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

    from scitbx import matrix
    B = self._xlucp.get_state()
    O = (B.transpose()).inverse()
    a, b, c = self._xlucp.get_model().get_unit_cell().parameters()[0:3]
    avec, bvec, cvec = [matrix.col(v) for v in O.transpose().as_list_of_lists()]

    # calculate d[B^T]/dp
    dB_dp = self._xlucp.get_ds_dp()
    dBT_dp = [dB.transpose() for dB in dB_dp]

    # calculate d[O]/dp
    dO_dp = [-O * dBT * O for dBT in dBT_dp]

    # get derivatives of angles wrt vectors
    def dangle(u, v):
      return [matrix.col(e) for e in angle_derivative_wrt_vectors(u,v)]
    dalpha_db, dalpha_dc = dangle(bvec, cvec)
    dbeta_da, dbeta_dc = dangle(avec, cvec)
    dgamma_da, dgamma_db = dangle(avec, bvec)

    da = []
    db = []
    dc = []
    daa = []
    dbb = []
    dcc = []

    # loop over parameters of the model
    for i, dO in enumerate(dO_dp):

      # extract derivatives of each unit cell vector wrt p
      dav_dp, dbv_dp, dcv_dp = dO.transpose().as_list_of_lists()
      dav_dp = matrix.col(dav_dp)
      dbv_dp = matrix.col(dbv_dp)
      dcv_dp = matrix.col(dcv_dp)

      # derivative of cell lengths wrt p
      da_dp = 1./a * avec.dot(dav_dp) if sel[0] else 0.0
      da.append(da_dp)
      db_dp = 1./b * bvec.dot(dbv_dp) if sel[1] else 0.0
      db.append(db_dp)
      dc_dp = 1./c * cvec.dot(dcv_dp) if sel[2] else 0.0
      dc.append(dc_dp)

      # derivative of cell angles wrt p
      daa_dp = dbv_dp.dot(dalpha_db) + dcv_dp.dot(dalpha_dc) if sel[3] else 0.0
      daa_dp *= RAD2DEG
      daa.append(daa_dp)
      dbb_dp = dav_dp.dot(dbeta_da) + dcv_dp.dot(dbeta_dc) if sel[4] else 0.0
      dbb_dp *= RAD2DEG
      dbb.append(dbb_dp)
      dcc_dp = dav_dp.dot(dgamma_da) + dbv_dp.dot(dgamma_db) if sel[5] else 0.0
      dcc_dp *= RAD2DEG
      dcc.append(dcc_dp)

    return da, db, dc, daa, dbb, dcc

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
    """For each residual, return the gradients dR/dp. Requires values to be
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

class GroupTie(object):
  """Base class for ties of multiple parameters together to a shared target
  value calculated on the group"""

  def __init__(self, parameter_list, weight):

    self._plist = parameter_list
    self._w = weight
    self._dRdp = None

  def values(self):

    raise NotImplementedError()

  def gradients(self):
    """Return dR/dp for each parameter."""
    return self._dRdp


class TiesToMean(GroupTie):
  """Tie a group of parameters together by similarity to their mean value"""

  # R = w Sum(p(i) - <p>)^2
  # dR/dp(i) = 2w(p(i) - <p>) d/dp(i)[ p(i) - <p>]
  # = 2w(p(i) - <p>) [1 - d<p>/dp(i)]
  # <p> = Sum(p(i))/n
  # d<p>/dp(i) = 1/n
  # dR/dp(i) = 2w(p(i) - <p>) [1 - 1/n]

  def __init__(self, *args, **kwargs):

    # set up base class
    super(TiestoMean, self).__init__(*args, **kwargs)

    # set scale factor for the gradient
    self._gradscal = 1.0 - 1.0/len(self._plist)

  def values(self):
    """Calculate and return squared residuals R for each parameter, cache
    gradients"""

    p = flex.double([p.value for p in self._plist])
    target = flex.mean(p)
    d = p - target
    wd = self._w * d
    self._dRdp = 2 * wd * self._gradscal

    return wd * d

class TiesToMedian(GroupTie):
  """Tie a group of parameters together by similarity to their median value"""

  def values(self):
    """Calculate and return squared residuals R for each parameter, cache
    gradients"""

    p = flex.double([p.value for p in self._plist])
    target = flex.median(p)
    d = p - target
    wd = self._w * d
    self._dRdp = 2 * wd # ignore dependence of median on the parameter values

    return wd * d

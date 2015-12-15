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
from dials.algorithms.refinement.refinement_helpers import AngleDerivativeWrtVectorElts
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

  def value(self, parameter_value, parameter_gradients):
    """Calculate and return weighted squared residual R, cache gradients"""

    d = parameter_value - self._target
    grad_coeff = 2. * self._w * d
    self._dRdp = [grad_coeff * g for g in parameter_gradients]

    return wd * d

  def gradient(self):
    """Return dR/dp"""
    return self._dRdp

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
    self._sigma = []
    for sig, grad in zip(sigma, grads):
      tst = [g > 1.e-10 for g in grad]
      if any(tst):
        self._sigma.append(sig)
      else:
        self._sigma.append(None)

    # For each non-zero sigma create a restraint between the relevant cell
    # parameter and its target value
    self._ties = []
    for t, s  in zip(self._target, self._sigma):
      if s is not None:
        self._ties.append(DerivedParameterTie(t, 1./s**2))
      else:
        self._ties.append(None)

    return

  def _calculate_uc_gradients(self, sel=[True]*6):
    '''Calculate gradients of the unit cell parameters with respect to '''

    from scitbx import matrix
    B = self._xlucp.get_state()
    O = (B.transpose()).inverse()
    a, b, c, aa, bb, cc = self._xlucp.get_model().get_unit_cell().parameters()
    aa *= DEG2RAD
    bb *= DEG2RAD
    cc *= DEG2RAD
    avec, bvec, cvec = self._xlucp.get_model().get_real_space_vectors()
    ua = avec.normalize()
    ub = bvec.normalize()
    uc = cvec.normalize()

    # calculate d[B^T]/dp
    dB_dp = self._xlucp.get_ds_dp()
    dBT_dp = [dB.transpose() for dB in dB_dp]

    # calculate d[O]/dp
    dO_dp = [-O * dBT * O for dBT in dBT_dp]

    # objects to get derivative of angles wrt vectors
    dalpha = AngleDerivativeWrtVectorElts(bvec, cvec)
    dbeta = AngleDerivativeWrtVectorElts(avec, cvec)
    dgamma = AngleDerivativeWrtVectorElts(avec, bvec)
    dalpha_db = dalpha.derivative_wrt_u()
    dalpha_dc = dalpha.derivative_wrt_v()
    dbeta_da = dbeta.derivative_wrt_u()
    dbeta_dc = dbeta.derivative_wrt_v()
    dgamma_da = dgamma.derivative_wrt_u()
    dgamma_db = dgamma.derivative_wrt_v()

    # XXX DEBUG compare with FD gradients
    fd_grads = self._check_fd_gradients()

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
      print "d[a]/dp{2} analytical: {0} FD: {1}".format(da_dp, fd_grads[i][0], i)

      db_dp = 1./b * bvec.dot(dbv_dp) if sel[1] else 0.0
      db.append(db_dp)
      print "d[b]/dp{2} analytical: {0} FD: {1}".format(db_dp, fd_grads[i][1], i)

      dc_dp = 1./c * cvec.dot(dcv_dp) if sel[2] else 0.0
      dc.append(dc_dp)
      print "d[c]/dp{2} analytical: {0} FD: {1}".format(dc_dp, fd_grads[i][2], i)

      # calculate orthogonal rate of change vectors
      if dav_dp.length() < 1e-10:
        ortho_dav_dp = matrix.col((0, 0, 0))
      else:
        v = avec.cross(dav_dp).normalize()
        u = v.cross(ua).normalize()
        ortho_dav_dp = dav_dp.dot(u) * u

      if dbv_dp.length() < 1e-10:
        ortho_dbv_dp = matrix.col((0, 0, 0))
      else:
        v = bvec.cross(dbv_dp).normalize()
        u = v.cross(ub).normalize()
        ortho_dbv_dp = dbv_dp.dot(u) * u

      if dcv_dp.length() < 1e-10:
        ortho_dcv_dp = matrix.col((0, 0, 0))
      else:
        v = cvec.cross(dcv_dp).normalize()
        u = v.cross(uc).normalize()
        ortho_dcv_dp = dcv_dp.dot(u) * u

      # derivative of cell angles wrt p
      daa_dp = ortho_dbv_dp.dot(dalpha_db) + ortho_dcv_dp.dot(dalpha_dc) \
        if sel[3] else 0.0
      daa.append(daa_dp)
      dbb_dp = ortho_dav_dp.dot(dbeta_da) + ortho_dcv_dp.dot(dbeta_dc) \
        if sel[4] else 0.0
      dbb.append(daa_dp)
      dcc_dp = ortho_dav_dp.dot(dgamma_da) + ortho_dbv_dp.dot(dgamma_db) \
        if sel[5] else 0.0
      dcc.append(dcc_dp)

      print "d[aa]/dp{2} analytical: {0} FD: {1}".format(RAD2DEG * daa_dp, fd_grads[i][3], i)
      print "d[bb]/dp{2} analytical: {0} FD: {1}".format(RAD2DEG * dbb_dp, fd_grads[i][4], i)
      print "d[cc]/dp{2} analytical: {0} FD: {1}".format(RAD2DEG * dcc_dp, fd_grads[i][5], i)

    return da, db, dc, daa, dbb, dcc

  def _check_fd_gradients(self):

    print "CHECKING GRADIENTS."
    from scitbx import matrix
    mp = self._xlucp
    p_vals = mp.get_param_vals()
    deltas = [1.e-7 for p in p_vals]
    assert len(deltas) == len(p_vals)
    fd_grad = []

    for i in range(len(deltas)):

      val = p_vals[i]

      p_vals[i] -= deltas[i] / 2.
      mp.set_param_vals(p_vals)
      rev_state = mp.get_model().get_unit_cell().parameters()

      p_vals[i] += deltas[i]
      mp.set_param_vals(p_vals)
      fwd_state = mp.get_model().get_unit_cell().parameters()

      fd_grad.append([(f - r) / deltas[i] for f,r in zip(fwd_state, rev_state)])

      p_vals[i] = val

    # return to the initial state
    mp.set_param_vals(p_vals)

    return fd_grad

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

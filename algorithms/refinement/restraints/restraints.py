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
DEG2RAD = pi/180.0
RAD2DEG = 180.0/pi

class SingleTie(object):
  """Tie of a single parameter to a value with a least squares restraint"""

  def __init__(self, parameter, target, weight):

    self._parameter = parameter
    self._target = target
    self._w = weight
    self._dRdp = None

  def value(self):
    """Calculate and return weighted squared residual R, cache gradient"""

    d = self._parameter.value - self._target
    wd = self._w * d
    self._dRdp = 2. * wd

    return wd * d

  def gradient(self):
    """Return dR/dp"""
    return self._dRdp

# FIXME SingleModelTie would be neater - express all restraints for parameters
# of a particular model in terms of the elements of the state of that model.
# calculate derivatives wrt parameters of the model using d[state]/dp elements.

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
    self._sigma = sigma

    assert len(self._target) == 6
    assert len(self._sigma) == 6

    # calculate gradients of cell parameters wrt model parameters. A gradient
    # of zero means that cell parameter is constrained and thus ignored in
    # restraints

    grads = self._calculate_uc_gradients()
    #self._sigma = [s if g > 1.e-12 else None for s,g in zip(self._sigma, grads)]

  def _calculate_uc_gradients(self):

    from scitbx import matrix
    B = self._xlucp.get_state()
    O = (B.transpose()).inverse()
    a, b, c, aa, bb, cc = self._xlucp.get_model().get_unit_cell().parameters()
    aa *= DEG2RAD
    bb *= DEG2RAD
    cc *= DEG2RAD
    avec, bvec, cvec = self._xlucp.get_model().get_real_space_vectors()

    # calculate d[B^T]/dp
    dB_dp = self._xlucp.get_ds_dp()
    dBT_dp = [dB.transpose() for dB in dB_dp]

    # calculate d[O]/dp
    dO_dp = [-O * dBT * O for dBT in dBT_dp]

    # XXX DEBUG compare with FD gradients
    fd_grads = self._check_fd_gradients()

    # loop over parameters of the model
    for i, dO in enumerate(dO_dp):

      # extract derivatives of each unit cell vector wrt p
      dav_dp, dbv_dp, dcv_dp = dO.transpose().as_list_of_lists()
      dav_dp = matrix.col(dav_dp)
      dbv_dp = matrix.col(dbv_dp)
      dcv_dp = matrix.col(dcv_dp)

      # derivative of cell params wrt p
      da_dp = 1./a * avec.dot(dav_dp)
      print "d[a]/dp{2} analytical: {0} FD: {1}".format(da_dp, fd_grads[i][0], i)

      db_dp = 1./b * bvec.dot(dbv_dp)
      print "d[b]/dp{2} analytical: {0} FD: {1}".format(db_dp, fd_grads[i][1], i)

      dc_dp = 1./c * cvec.dot(dcv_dp)
      print "d[c]/dp{2} analytical: {0} FD: {1}".format(dc_dp, fd_grads[i][2], i)

      z = bvec.dot(cvec) / (b * c)
      daa_dp = bvec.dot(cvec) * (db_dp * c + b * dc_dp) - b * c * (dbv_dp.dot(cvec) + bvec.dot(dcv_dp))
      daa_dp /= (b * b * c * c)
      daa_dp *= -RAD2DEG / (sqrt(1 - z**2))
      #daa_dp = -1. * cos(aa) * (db_dp * c + b * dc_dp)
      #daa_dp += dbv_dp.dot(cvec) + bvec.dot(dcv_dp)
      #daa_dp *= RAD2DEG / (b*c*sin(aa))
      print "d[alpha]/dp{2} analytical: {0} FD: {1}".format(daa_dp, fd_grads[i][3], i)

    #from dials.util.command_line import interactive_console; interactive_console()
    #1/0

  def _check_fd_gradients(self):

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

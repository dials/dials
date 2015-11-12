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

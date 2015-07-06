#!/usr/bin/env python
#
#  __init__.py
#
#  Copyright (C) 2015 Diamond Light Source and STFC Rutherford Appleton
#                     Laboratory, UK.
#
#  Author: David Waterman
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division
from logging import info, debug

class CentroidOutlier(object):
  """Base class for centroid outlier detection algorithms"""

  def __init__(self, cols=["x_resid", "y_resid", "phi_resid"],
               min_num_obs=20):

    # minimum number of observations per panel below which all reflections will
    # be marked as potential outliers
    self._min_num_obs = min_num_obs

    # the number of rejections
    self.nreject = 0

    return

  def _detect_outliers(cols):
    """Perform outlier detection using the input cols and return a flex.bool
    indicating which rows in the cols are considered outlying. cols should be
    a list of flex.doubles of equal lengths"""

    # to be implemented by derived classes
    raise NotImplementedError()

  def __call__(self, reflections):
    """Identify outliers in the input and set the centroid_outlier flag.
    Return True if any outliers were detected, otherwise False"""

    # check the columns are present
    keys = reflections.keys()
    for col in cols: assert col in keys
    self._cols = cols

    sel = reflections.get_flags(reflections.flags.used_in_refinement)
    matches = reflections.select(sel)
    all_imatches = sel.iselection()

    max_panel = flex.max(matches['panel'])

    nreject = 0
    for pnl in xrange(max_panel + 1):

      # selection of reflections that are on this panel
      pnl_sel = matches['panel'] == pnl

      # their indices back in the original reflection table
      pnl_isel = all_imatches.select(pnl_sel)

      if len(pnl_isel) >= self._min_num_obs:

        # get the subset of data on this panel as a list of columns
        cols = [matches[col].select(pnl_sel) for col in self._cols]

        # determine the position of outliers on this panel's sub-dataset
        outliers = self._detect_outliers(cols)

        # get positions of outliers from the original matches
        ioutliers = pnl_isel.select(outliers)

      else:
        msg = "Only {0} reflections on panel {1}. ".format(len(pnl_isel), pnl)
        debug(msg + "All of these flagged as possible outliers.")
        ioutliers = pnl_isel

      # set those reflections as outliers in the original reflection table
      reflections.set_flags(ioutliers,
        reflections.flags.centroid_outlier)

      self.nreject += len(ioutliers)

    if self.nreject == 0: return False

    info("{0} reflections have been flagged as outliers".format(nreject))

    return True


#!/usr/bin/env python
#
#  tukey.py
#
#  Copyright (C) 2015 Diamond Light Source and STFC Rutherford Appleton
#                     Laboratory, UK.
#
#  Author: David Waterman
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division
from dials.algorithms.refinement.outlier_detection import CentroidOutlier
from dials.array_family import flex

class Tukey(CentroidOutlier):
  """Implementation of the CentroidOutlier class using Tukey's rule of thumb.
  That is values more than iqr_multiplier times the interquartile range from
  the quartiles are designed outliers. When x=1.5, this is Tukey's rule."""

  def __init__(self, cols=["x_resid", "y_resid", "phi_resid"],
               min_num_obs=20,
               separate_experiments=True,
               separate_panels=True,
               iqr_multiplier=1.5):

    CentroidOutlier.__init__(self,
      cols=["x_resid", "y_resid", "phi_resid"],
      min_num_obs=20,
      separate_experiments=True,
      separate_panels=True)

    self._iqr_multiplier = iqr_multiplier

    return

  def _detect_outliers(self, cols):

    from scitbx.math import five_number_summary

    outliers = flex.bool(len(cols[0]), False)
    for col in cols:
      min_x, q1_x, med_x, q3_x, max_x = five_number_summary(col)
      iqr_x = q3_x - q1_x
      cut_x = self._iqr_multiplier * iqr_x
      outliers.set_selected(col > q3_x + cut_x, True)
      outliers.set_selected(col < q1_x - cut_x, True)

    return outliers

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

from __future__ import absolute_import, division

from dials.algorithms.refinement.outlier_detection import CentroidOutlier
from dials.array_family import flex

class SauterPoon(CentroidOutlier):
  """Implementation of the CentroidOutlier class using the algorithm of
  Sauter & Poon (2010) (https://doi.org/10.1107/S0021889810010782)."""

  def __init__(self, cols=None,
               min_num_obs=20,
               separate_experiments=True,
               separate_panels=True,
               block_width=None,
               px_sz=(1,1),
               verbose=False,
               pdf=None):

    # here the column names are fixed by the algorithm, so what's passed in is
    # ignored.
    CentroidOutlier.__init__(self,
      cols=["miller_index", "xyzobs.px.value", "xyzcal.px"],
      min_num_obs=min_num_obs,
      separate_experiments=separate_experiments,
      separate_panels=separate_panels,
      block_width=block_width)

    self._px_sz = px_sz
    self._verbose = verbose
    self._pdf = pdf

    return

  def _detect_outliers(self, cols):

    # cols is guaranteed to be a list of three flex arrays, containing miller
    # indices, observed pixel coordinates and calculated pixel coordinates.
    # Copy the data into matches
    class match: pass
    matches = []
    for hkl in cols[0]:
      m = match()
      m.miller_index = hkl
      matches.append(m)

    for obs, m in zip(cols[1], matches):
      m.x_obs = obs[0]*self._px_sz[0]
      m.y_obs = obs[1]*self._px_sz[1]

    for calc, m in zip(cols[2], matches):
      m.x_calc = calc[0]*self._px_sz[0]
      m.y_calc = calc[1]*self._px_sz[1]

    from rstbx.phil.phil_preferences import indexing_api_defs
    import iotbx.phil
    hardcoded_phil = iotbx.phil.parse(
    input_string=indexing_api_defs).extract()

    # set params into the hardcoded_phil
    hardcoded_phil.indexing.outlier_detection.verbose = self._verbose
    hardcoded_phil.indexing.outlier_detection.pdf = self._pdf

    from rstbx.indexing_api.outlier_procedure import OutlierPlotPDF
    if self._pdf is not None:
    ## new code for outlier rejection inline here
      hardcoded_phil.__inject__("writer",OutlierPlotPDF(hardcoded_phil.indexing.outlier_detection.pdf))

    # execute Sauter and Poon (2010) algorithm
    from rstbx.indexing_api import outlier_detection
    hardcoded_phil = hardcoded_phil
    od = outlier_detection.find_outliers_from_matches(
      matches,
      verbose=self._verbose,
      horizon_phil=hardcoded_phil)

    # flex.bool of the inliers
    outliers = ~od.get_cache_status()

    return outliers

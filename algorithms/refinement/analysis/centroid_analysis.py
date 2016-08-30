#!/usr/bin/env cctbx.python
#
#  Copyright (C) (2016) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman.
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#
"""Analysis of centroid residuals for determining suitable refinement and
outlier rejection parameters automatically"""

from __future__ import division
from math import pi, floor, ceil
from dials.array_family import flex
from periodogram import Periodogram

RAD2DEG = 180./pi

class CentroidAnalyser(object):

  def __init__(self, reflections, av_callback=flex.mean):

    self._av_callback = av_callback

    # Remove invalid/bad reflections
    reflections = reflections.select(reflections.get_flags(
      reflections.flags.indexed))
    reflections = reflections.select(~(reflections['miller_index'] == (0,0,0)))

    # Trim extrema of the scan as they have the worst phi residuals
    phi_obs = reflections['xyzobs.mm.value'].parts()[2]
    phi_max, phi_min = flex.max(phi_obs), flex.min(phi_obs)
    sel = (phi_obs > phi_min + 0.005) & (phi_obs < phi_max - 0.005)
    reflections = reflections.select(sel)

    # FIXME - better way to recognise non-predictions. Can't rely on flags
    # in e.g. indexed.pickle I think.
    x, y, z = reflections['xyzcal.mm'].parts()
    sel = (x == 0) & (y == 0)
    reflections = reflections.select(~sel)
    self._nexp = flex.max(reflections['id'])

    # Ensure required keys are present
    if not all([k in reflections for k in ['x_resid', 'y_resid', 'phi_resid']]):
      x_obs, y_obs, phi_obs = reflections['xyzobs.mm.value'].parts()
      x_cal, y_cal, phi_cal = reflections['xyzcal.mm'].parts()

      # do not wrap around multiples of 2*pi; keep the full rotation
      # from zero to differentiate repeat observations.
      from math import pi
      TWO_PI = 2.0 * pi
      resid = phi_cal - (flex.fmod_positive(phi_obs, TWO_PI))
      # ensure this is the smaller of two possibilities
      resid = flex.fmod_positive((resid + pi), TWO_PI) - pi
      phi_cal = phi_obs + resid
      reflections['x_resid'] = x_cal - x_obs
      reflections['y_resid'] = y_cal - y_obs
      reflections['phi_resid'] = phi_cal - phi_obs

    self._reflections = reflections

  def __call__(self):
    """Perform power spectrum analysis and return the results as a list
    of dictionaries (one for each experiment)"""

    results = []
    for iexp in range(self._nexp + 1):
      reflections = self._reflections.select(self._reflections['id'] == iexp)
      x_resid = reflections['x_resid']
      y_resid = reflections['y_resid']
      phi_resid = reflections['phi_resid']
      phi_obs_deg = reflections['xyzobs.mm.value'].parts()[2] * RAD2DEG

      # Calculate average residuals in equal-width blocks
      phi_range = flex.min(phi_obs_deg), flex.max(phi_obs_deg)
      phi_width = phi_range[1] - phi_range[0]
      ideal_block_size = 1.0
      while True:
        nblocks = int(phi_width // ideal_block_size)
        block_size = phi_width / nblocks
        xr_per_blk = flex.double()
        yr_per_blk = flex.double()
        pr_per_blk = flex.double()
        nr = flex.int()
        for i in range(nblocks - 1):
          blk_start = phi_range[0] + i * block_size
          blk_end = blk_start + block_size
          sel = (phi_obs_deg >= blk_start) & (phi_obs_deg < blk_end)
          nr.append(sel.count(True))
          xr_per_blk.append(self._av_callback(x_resid.select(sel)))
          yr_per_blk.append(self._av_callback(y_resid.select(sel)))
          pr_per_blk.append(self._av_callback(phi_resid.select(sel)))
        # include max phi in the final block
        blk_start = phi_range[0] + (nblocks - 1) * block_size
        blk_end = phi_range[1]
        sel = (phi_obs_deg >= blk_start) & (phi_obs_deg <= blk_end)
        nr.append(sel.count(True))
        xr_per_blk.append(self._av_callback(x_resid.select(sel)))
        yr_per_blk.append(self._av_callback(y_resid.select(sel)))
        pr_per_blk.append(self._av_callback(phi_resid.select(sel)))

        # Break if there are enough reflections, otherwise increase block size
        min_nr = flex.min(nr)
        if min_nr > 50: break
        fac = 50 / min_nr
        ideal_block_size *= fac

      results_this_exp = {'block_size':block_size}

      # Perform power spectrum analysis on the residuals
      px = Periodogram(xr_per_blk)
      # FIXME here extract information from the power spectrum

      # collect results
      results.append(results_this_exp)

      print "block_size", block_size
      print "#refs per block"
      print list(nr)

      import matplotlib.pyplot as plt
      plt.plot(xr_per_blk)
      plt.ylabel('some numbers')
      plt.show()

    return results

if __name__ == "__main__":

  import sys
  ref = sys.argv[1]

  from dials.array_family import flex
  refs = flex.reflection_table.from_pickle(ref)

  ca = CentroidAnalyser(refs)
  ca()

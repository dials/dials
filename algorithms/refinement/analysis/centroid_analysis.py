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

    self._spectral_analysis = False
    self._av_callback = av_callback

    # Remove invalid reflections
    reflections = reflections.select(~(reflections['miller_index'] == (0,0,0)))

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

    # start populating results dictionary
    self._results = []
    for iexp in range(self._nexp + 1):
      ref_this_exp = reflections.select(reflections['id'] == iexp)
      if len(ref_this_exp) == 0:
        self._results.append({})
        continue
      x_resid = ref_this_exp['x_resid']
      y_resid = ref_this_exp['y_resid']
      phi_resid = ref_this_exp['phi_resid']
      phi_obs_deg = ref_this_exp['xyzobs.mm.value'].parts()[2] * RAD2DEG

      # Calculate average residuals in equal-width blocks
      phi_range = flex.min(phi_obs_deg), flex.max(phi_obs_deg)
      phi_width = phi_range[1] - phi_range[0]
      ideal_block_size = 1.0
      while True:
        nblocks = max(int(phi_width // ideal_block_size), 1)
        block_size = phi_width / nblocks
        xr_per_blk = flex.double()
        yr_per_blk = flex.double()
        pr_per_blk = flex.double()
        nr = flex.int()
        for i in range(nblocks - 1):
          blk_start = phi_range[0] + i * block_size
          blk_end = blk_start + block_size
          sel = (phi_obs_deg >= blk_start) & (phi_obs_deg < blk_end)
          nref_in_block = sel.count(True)
          nr.append(nref_in_block)
          if nref_in_block > 0:
            xr_per_blk.append(self._av_callback(x_resid.select(sel)))
            yr_per_blk.append(self._av_callback(y_resid.select(sel)))
            pr_per_blk.append(self._av_callback(phi_resid.select(sel)))
          else:
            xr_per_blk.append(0.0)
            yr_per_blk.append(0.0)
            pr_per_blk.append(0.0)
        # include max phi in the final block
        blk_start = phi_range[0] + (nblocks - 1) * block_size
        blk_end = phi_range[1]
        sel = (phi_obs_deg >= blk_start) & (phi_obs_deg <= blk_end)
        nref_in_block = sel.count(True)
        nr.append(nref_in_block)
        if nref_in_block > 0:
          xr_per_blk.append(self._av_callback(x_resid.select(sel)))
          yr_per_blk.append(self._av_callback(y_resid.select(sel)))
          pr_per_blk.append(self._av_callback(phi_resid.select(sel)))
        else:
            xr_per_blk.append(0.0)
            yr_per_blk.append(0.0)
            pr_per_blk.append(0.0)

        # Break if there are enough reflections, otherwise increase block size,
        # unless only one block remains
        if nblocks == 1: break
        min_nr = flex.min(nr)
        if min_nr > 50: break
        if min_nr < 5:
          fac = 2
        else:
          fac = 50 / min_nr
        ideal_block_size *= fac

      results_this_exp = {'block_size':block_size,
                          'nref_per_block':nr,
                          'av_x_resid_per_block':xr_per_blk,
                          'av_y_resid_per_block':yr_per_blk,
                          'av_phi_resid_per_block':pr_per_blk,}

      # collect results
      self._results.append(results_this_exp)

  def __call__(self, do_spectral_analysis=True):
    """Perform analysis and return the results as a list of dictionaries (one
    for each experiment)"""

    if do_spectral_analysis:
      if self._spectral_analysis: return self._results

      # Perform power spectrum analysis on the residuals
      for exp_data in self._results:
        px = Periodogram(exp_data['av_x_resid_per_block'])
        exp_data['x_periodogram'] = px
        py = Periodogram(exp_data['av_y_resid_per_block'])
        exp_data['y_periodogram'] = py
        pz = Periodogram(exp_data['av_phi_resid_per_block'])
        exp_data['phi_periodogram'] = pz

        # FIXME here extract further information from the power spectrum

    return self._results

if __name__ == "__main__":

  import sys
  import matplotlib.pyplot as plt
  from dials.array_family import flex

  ref = sys.argv[1]
  refs = flex.reflection_table.from_pickle(ref)

  ca = CentroidAnalyser(refs)
  results = ca()

  for e in results:

    plt.plot(e['av_x_resid_per_block'])
    plt.ylabel('x residuals per block')
    plt.show()
    e['x_periodogram'].plot()

    plt.plot(e['av_y_resid_per_block'])
    plt.ylabel('y residuals per block')
    plt.show()
    e['y_periodogram'].plot()

    plt.plot(e['av_phi_resid_per_block'])
    plt.ylabel('phi residuals per block')
    plt.show()
    e['phi_periodogram'].plot()

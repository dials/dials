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

    # flags to indicate at what level the analysis has been performed
    self._average_residuals = False
    self._spectral_analysis = False

    self._av_callback = av_callback

    # Remove invalid reflections
    reflections = reflections.select(~(reflections['miller_index'] == (0,0,0)))

    # FIXME - better way to recognise non-predictions. Can't rely on flags
    # in e.g. indexed.pickle I think.
    x, y, z = reflections['xyzcal.mm'].parts()
    sel = (x == 0) & (y == 0)
    reflections = reflections.select(~sel)
    self._nexp = flex.max(reflections['id']) + 1

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

    # create empty results list
    self._results = []

    # first, just determine a suitable blocksize for analysis
    for iexp in range(self._nexp):
      ref_this_exp = reflections.select(reflections['id'] == iexp)
      if len(ref_this_exp) == 0:
        # can't do anything, just keep an empty dictionary
        self._results.append({})
        continue
      phi_obs_deg = ref_this_exp['xyzobs.mm.value'].parts()[2] * RAD2DEG
      phi_range = flex.min(phi_obs_deg), flex.max(phi_obs_deg)
      phi_width = phi_range[1] - phi_range[0]
      ideal_block_size = 1.0
      old_nblocks = 0
      while True:
        nblocks = int(phi_width // ideal_block_size)
        if nblocks == old_nblocks: nblocks -= 1
        nblocks = max(nblocks, 1)
        block_size = phi_width / nblocks
        nr = flex.int()
        for i in range(nblocks - 1):
          blk_start = phi_range[0] + i * block_size
          blk_end = blk_start + block_size
          sel = (phi_obs_deg >= blk_start) & (phi_obs_deg < blk_end)
          nref_in_block = sel.count(True)
          nr.append(nref_in_block)
        # include max phi in the final block
        blk_start = phi_range[0] + (nblocks - 1) * block_size
        blk_end = phi_range[1]
        sel = (phi_obs_deg >= blk_start) & (phi_obs_deg <= blk_end)
        nref_in_block = sel.count(True)
        nr.append(nref_in_block)
        # Break if there are enough reflections, otherwise increase block size,
        # unless only one block remains
        if nblocks == 1: break
        min_nr = flex.min(nr)
        if min_nr >= 50: break
        if min_nr < 5:
          fac = 2
        else:
          fac = 50 / min_nr
        ideal_block_size *= fac
        old_nblocks = nblocks

      # collect the basic data for this experiment
      self._results.append({'block_size':block_size,
                            'nref_per_block':nr,
                            'nblocks':nblocks,
                            'phi_range':phi_range})

    # keep reflections for analysis
    self._reflections = reflections

  def __call__(self, calc_average_residuals=True,
                     calc_periodograms=True, spans=[4,4]):
    """Perform analysis and return the results as a list of dictionaries (one
    for each experiment)"""

    # if not doing further analysis, return the basic data
    if not calc_average_residuals and not calc_periodograms:
      return self._results

    # if we don't have average residuals already, calculate them
    if not self._average_residuals:
      for iexp in range(self._nexp):
        results_this_exp = self._results[iexp]
        block_size = results_this_exp.get('block_size')
        if block_size is None: continue
        phi_range = results_this_exp['phi_range']
        phi_width = phi_range[1] - phi_range[0]
        nblocks = results_this_exp['nblocks']
        block_size = phi_width / nblocks
        ref_this_exp = self._reflections.select(self._reflections['id'] == iexp)
        x_resid = ref_this_exp['x_resid']
        y_resid = ref_this_exp['y_resid']
        phi_resid = ref_this_exp['phi_resid']
        phi_obs_deg = ref_this_exp['xyzobs.mm.value'].parts()[2] * RAD2DEG
        xr_per_blk = flex.double()
        yr_per_blk = flex.double()
        pr_per_blk = flex.double()
        for i in range(nblocks - 1):
          blk_start = phi_range[0] + i * block_size
          blk_end = blk_start + block_size
          sel = (phi_obs_deg >= blk_start) & (phi_obs_deg < blk_end)
          xr_per_blk.append(self._av_callback(x_resid.select(sel)))
          yr_per_blk.append(self._av_callback(y_resid.select(sel)))
          pr_per_blk.append(self._av_callback(phi_resid.select(sel)))
        # include max phi in the final block
        blk_start = phi_range[0] + (nblocks - 1) * block_size
        blk_end = phi_range[1]
        sel = (phi_obs_deg >= blk_start) & (phi_obs_deg <= blk_end)
        xr_per_blk.append(self._av_callback(x_resid.select(sel)))
        yr_per_blk.append(self._av_callback(y_resid.select(sel)))
        pr_per_blk.append(self._av_callback(phi_resid.select(sel)))

        results_this_exp['av_x_resid_per_block'] = xr_per_blk
        results_this_exp['av_y_resid_per_block'] = yr_per_blk
        results_this_exp['av_phi_resid_per_block'] = pr_per_blk
      self._average_residuals = True

    # Perform power spectrum analysis on the residuals
    if calc_periodograms:
      if self._spectral_analysis: return self._results

      for exp_data in self._results:
        px = Periodogram(exp_data['av_x_resid_per_block'], spans=spans)
        exp_data['x_periodogram'] = px
        py = Periodogram(exp_data['av_y_resid_per_block'], spans=spans)
        exp_data['y_periodogram'] = py
        pz = Periodogram(exp_data['av_phi_resid_per_block'], spans=spans)
        exp_data['phi_periodogram'] = pz
      self._spectral_analysis = True

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

    nblocks = e['nblocks']
    block_size = e['block_size']
    phistart = e['phi_range'][0]
    block_centres = block_size * flex.double_range(nblocks) + phistart + block_size/2.0
    plt.plot(block_centres, e['av_x_resid_per_block'])
    plt.xlabel('phi (degrees)')
    plt.ylabel('x residuals per block')
    plt.show()
    e['x_periodogram'].plot(sample_interval=block_size)

    plt.plot(block_centres, e['av_y_resid_per_block'])
    plt.xlabel('phi (degrees)')
    plt.ylabel('y residuals per block')
    plt.show()
    e['y_periodogram'].plot(sample_interval=block_size)

    plt.plot(block_centres, e['av_phi_resid_per_block'])
    plt.xlabel('phi (degrees)')
    plt.ylabel('phi residuals per block')
    plt.show()
    e['phi_periodogram'].plot(sample_interval=block_size)

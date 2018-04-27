from __future__ import absolute_import, division
from __future__ import print_function
#!/usr/bin/env python
# -*- mode: python; coding: utf-8; indent-tabs-mode: nil; python-indent: 2 -*-
#
# dials.algorithms.indexing.refinement.py
#
#  Copyright (C) 2014 Diamond Light Source
#
#  Author: Richard Gildea
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from cctbx.array_family import flex
import logging
import math

logger = logging.getLogger(__name__)

def refine(params, reflections, experiments,
           verbosity=0, debug_plots=False):
  detector = experiments.detectors()[0]
  if params.refinement.parameterisation.scan_varying:
    logger.warn(
      'scan_varying=True not supported in indexing: setting scan_varying=False')
    params.refinement.parameterisation.scan_varying = False

  from dials.algorithms.refinement import RefinerFactory
  refiner = RefinerFactory.from_parameters_data_experiments(
    params, reflections, experiments, verbosity=verbosity)

  outliers = None
  refined = refiner.run()
  if debug_plots:
    debug_plot_residuals(refiner)
  return refiner, refined, outliers


def debug_plot_residuals(refiner, inlier_sel=None):
  from matplotlib import pyplot
  matches = refiner.get_matches()
  x_residuals = matches['x_resid']
  y_residuals = matches['y_resid']
  phi_residuals = matches['phi_resid']
  residuals = flex.vec3_double(x_residuals, y_residuals, phi_residuals)
  frame_obs = matches['xyzobs.px.value'].parts()[2]
  phi_obs = matches['xyzobs.mm.value'].parts()[2]
  panel_ids = matches['panel']
  crystal_ids = matches['id']
  if inlier_sel is None:
    inlier_sel = flex.bool(len(residuals), True)
  print(inlier_sel.size(), panel_ids.size())
  for i_crystal in range(flex.max(crystal_ids)+1):
    crystal_sel = (crystal_ids == i_crystal)
    for i_panel in range(len(refiner.get_experiments().detectors()[0])):
      panel_sel = (panel_ids == i_panel)

      pyplot.axhline(0, color='grey')
      pyplot.axvline(0, color='grey')

      pyplot.scatter(
        x_residuals.select(inlier_sel & panel_sel & crystal_sel).as_numpy_array(),
        y_residuals.select(inlier_sel & panel_sel & crystal_sel).as_numpy_array(),
        c='b', alpha=0.5)
      pyplot.scatter(
        x_residuals.select((~inlier_sel) & panel_sel & crystal_sel).as_numpy_array(),
        y_residuals.select((~inlier_sel) & panel_sel & crystal_sel).as_numpy_array(),
        c='r', alpha=0.5)
      pyplot.axes().set_aspect('equal')
      pyplot.show()

  min_frame = int(math.floor(flex.min(frame_obs)))
  max_frame = int(math.ceil(flex.max(frame_obs)))
  mean_residuals_x = []
  mean_residuals_y = []
  mean_residuals_phi = []
  frame = []
  phi_obs_deg = (180/math.pi) * phi_obs
  phi = []

  for i_phi in range(int(math.floor(flex.min(phi_obs_deg))),
                 int(math.ceil(flex.max(phi_obs_deg)))):
    sel = (phi_obs_deg >= i_phi) & (phi_obs_deg < (i_phi+1))
    if sel.count(True) == 0:
      continue
    mean_residuals_x.append(flex.mean(x_residuals.select(sel)))
    mean_residuals_y.append(flex.mean(y_residuals.select(sel)))
    mean_residuals_phi.append(flex.mean(phi_residuals.select(sel)))
    phi.append(i_phi)

  fig = pyplot.figure()
  ax = fig.add_subplot(311)
  pyplot.axhline(0, color='grey')
  ax.scatter(phi, mean_residuals_x)
  ax.set_xlabel('phi (deg)')
  ax.set_ylabel('mean residual_x')
  ax = fig.add_subplot(312)
  pyplot.axhline(0, color='grey')
  ax.scatter(phi, mean_residuals_y)
  ax.set_xlabel('phi (deg)')
  ax.set_ylabel('mean residual_y')
  ax = fig.add_subplot(313)
  pyplot.axhline(0, color='grey')
  ax.scatter(phi, mean_residuals_phi)
  ax.set_xlabel('phi (deg)')
  ax.set_ylabel('mean residual_phi')
  pyplot.show()

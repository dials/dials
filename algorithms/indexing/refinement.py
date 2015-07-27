from __future__ import division
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

from libtbx.utils import plural_s
from cctbx.array_family import flex
import math


def refine(params, reflections, experiments, maximum_spot_error=None,
           maximum_phi_error=None,
           verbosity=0, debug_plots=False):
  detector = experiments.detectors()[0]

  from dials.algorithms.refinement import RefinerFactory
  refiner = RefinerFactory.from_parameters_data_experiments(
    params, reflections, experiments, verbosity=verbosity)

  outliers = None

  if maximum_spot_error is not None or maximum_phi_error is not None:
    matches = refiner.get_matches()
    x_residuals = matches['x_resid']
    y_residuals = matches['y_resid']
    phi_residuals = matches['phi_resid']
    residuals = flex.vec3_double(x_residuals, y_residuals, phi_residuals)
    frame_obs = matches['xyzobs.px.value'].parts()[2]
    panel_ids = matches['panel']
    crystal_ids = matches['id']
    match_iobs = matches['iobs']
    mm_residual_norms = flex.sqrt(
      flex.pow2(x_residuals) + flex.pow2(y_residuals))
    inlier_sel = flex.bool(len(matches), True)
    if maximum_spot_error is not None:
      # hard cutoff, but this is essentially what XDS does by default
      # assumes pixel size is same for all panels and same in x and y
      inlier_sel &= mm_residual_norms < (
        maximum_spot_error * detector[0].get_pixel_size()[0])
    if maximum_phi_error is not None:
      inlier_sel &= flex.abs(phi_residuals) < (math.pi * maximum_phi_error/180)

    print "Rejecting %i outlier%s" %plural_s(inlier_sel.count(False))
    if debug_plots:
      debug_plot_residuals(refiner, inlier_sel=inlier_sel)

    # XXX Hack to do the outlier rejection without instatiating a new refiner
    # XXX TODO move this outlier rejection into the refinement proper
    refiner._refman._reflections = refiner._refman._reflections.select(inlier_sel)

    # sort the inliers so they are in the same order
    perm = flex.sort_permutation(match_iobs)
    outliers = (~inlier_sel.select(perm)).iselection()

    #reflections_for_refinement = reflections_for_refinement.select(
      #refiner.selection_used_for_refinement()).select(inlier_sel.select(perm))
    #refiner = RefinerFactory.from_parameters_data_experiments(
      #params, reflections_for_refinement, experiments,
      #verbosity=verbosity)

  matches = refiner.get_matches()
  crystal_ids = matches['id']
  # DGW commented out as reflections.minimum_number_of_reflections no longer exists
  #for i_cryst in range(flex.max(crystal_ids) + 1):
  #  if (crystal_ids == i_cryst).count(True) < params.refinement.reflections.minimum_number_of_reflections:
  #    raise RuntimeError("Insufficient matches for crystal %i" %(i_cryst+1))
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
  print inlier_sel.size(), panel_ids.size()
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
      #r = maximum_spot_error * self.detector[0].get_pixel_size()
      #pyplot.Circle((r, r), 0.5, color='b', fill=False)
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

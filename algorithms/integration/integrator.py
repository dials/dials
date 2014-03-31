#!/usr/bin/env python
#
# dials.algorithms.integration.integrator.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division

class Integrator(object):
  ''' The integrator base class. '''

  def __init__(self, n_sigma, n_blocks, filter_by_zeta):
    ''' Initialise the integrator base class.

    Params:

    '''
    self.n_sigma = n_sigma
    self.n_blocks = n_blocks
    self.filter_by_zeta = filter_by_zeta

  def __call__(self, experiments, reference=None, extracted=None):
    ''' Call to integrate.

    Params:
        sweep The sweep to process
        crystal The crystal to process
        reflections The reflection list
        reference The reference profiles

    Returns:
        A reflection list

    '''
    from dials.algorithms.shoebox import ReflectionBlockExtractor
    from dials.array_family import flex
    from dials.algorithms.shoebox import MaskCode

    assert(len(experiments) == 1)

    # Predict a load of reflections
    if extracted == None:
      predicted = flex.reflection_table.from_predictions(experiments[0])
      predicted['id'] = flex.size_t(len(predicted), 0)
    else:
      predicted = None

    from dials.framework.registry import Registry
    registry = Registry()
    params = registry.params()
    if params.integration.shoebox.sigma_b is None or params.integration.shoebox.sigma_m is None:
      assert(reference is not None)
      from dials.algorithms.profile_model.profile_model import ProfileModel
      from math import pi
      profile_model = ProfileModel(experiments[0], reference)
      params.integration.shoebox.sigma_b = profile_model.sigma_b() * 180.0 / pi
      params.integration.shoebox.sigma_m = profile_model.sigma_m() * 180.0 / pi
      print 'Sigma B: %f' % params.integration.shoebox.sigma_b
      print 'Sigma M: %f' % params.integration.shoebox.sigma_m

    # Get the extractor
    extract = ReflectionBlockExtractor(experiments[0], predicted,
      self.n_sigma, self.n_blocks, self.filter_by_zeta, extracted)

    # Loop through all the blocks
    result = flex.reflection_table()
    print ''
    for reflections in extract:
      from dials.algorithms import filtering
      from dials.util.command_line import Command
      # Set all reflections which overlap bad pixels to zero
      Command.start('Filtering reflections by detector mask')
      array_range = experiments[0].scan.get_array_range()
      mask = filtering.by_detector_mask(
        reflections['bbox'], experiments[0].imageset[0] >= 0, array_range)
      reflections.del_selected(mask != True)
      Command.end('Filtered {0} reflections by detector mask'.format(
        len(reflections)))

      # Filter the reflections by zeta
      if self.filter_by_zeta > 0:
        Command.start('Filtering reflections by zeta >= {0}'.format(
            self.filter_by_zeta))
        mask = filtering.by_zeta(experiments[0].goniometer, experiments[0].beam,
            reflections['s1'], self.filter_by_zeta)
        reflections.del_selected(mask != True)
        Command.end('Filtered {0} reflections by zeta >= {1}'.format(
          len(reflections), self.filter_by_zeta))

      if reference:
        from dials.algorithms.peak_finding.spot_matcher import SpotMatcher
        match = SpotMatcher(max_separation=1)
        sind, pind = match(reference, reflections)
        reflections.set_flags(pind, reflections.flags.reference_spot)

      reflections.integrate(experiments[0])

      # Compute number of foreground and background pixels
      bg_code = MaskCode.Valid | MaskCode.BackgroundUsed
      fg_code = MaskCode.Valid | MaskCode.Foreground
      n_bg = reflections['shoebox'].count_mask_values(bg_code)
      n_fg = reflections['shoebox'].count_mask_values(fg_code)
      reflections['n_background'] = n_bg
      reflections['n_foreground'] = n_fg

      # Delete the profiles
      del reflections['shoebox']
      del reflections['rs_shoebox']
      result.extend(reflections)
      print ''

    # Return the reflections
    result.sort('miller_index')
    return result

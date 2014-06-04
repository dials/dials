#
# integrator.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division


class ReflectionBlockIntegrator(object):
  ''' A class to perform the integration. '''

  def __init__(self, params, experiments, extractor=None):
    ''' Initialise the integrator. '''
    from math import pi
    from dials.algorithms import shoebox

    # Ensure we have 1 experiment at the moment
    assert(len(experiments) == 1)
    assert(extractor is not None)

    # Save the parameters
    self.params = params
    self.experiments = experiments
    self.extractor = extractor

    # Create the shoebox masker
    n_sigma = params.integration.shoebox.n_sigma
    sigma_b = params.integration.shoebox.sigma_b
    sigma_m = params.integration.shoebox.sigma_m
    assert(n_sigma > 0)
    assert(sigma_b > 0)
    assert(sigma_m > 0)
    delta_b = n_sigma * sigma_b * pi / 180.0
    delta_m = n_sigma * sigma_m * pi / 180.0
    self._mask_profiles = shoebox.Masker(experiments[0], delta_b, delta_m)

  def integrate(self):
    ''' Integrate all the reflections. '''
    from dials.array_family import flex
    from dials.algorithms.shoebox import MaskCode
    result = flex.reflection_table()
    for indices, reflections in self.extractor:
      self._mask_profiles(reflections, None)
      reflections.integrate(self.experiments[0])
      bg_code = MaskCode.Valid | MaskCode.BackgroundUsed
      fg_code = MaskCode.Valid | MaskCode.Foreground
      n_bg = reflections['shoebox'].count_mask_values(bg_code)
      n_fg = reflections['shoebox'].count_mask_values(fg_code)
      reflections['n_background'] = n_bg
      reflections['n_foreground'] = n_fg
      del reflections['shoebox']
      del reflections['rs_shoebox']
      result.extend(reflections)
    assert(len(result) > 0)
    result.sort('miller_index')
    return result


class Integrator(object):
  ''' Integrate reflections '''

  def __init__(self, params, exlist, reference=None,
               predicted=None, shoeboxes=None):
    '''Initialise the script.'''

    # Load the extractor based on the input
    if shoeboxes is not None:
      extractor = self._load_extractor(shoeboxes, params, exlist)
    else:
      if reference:
        self._compute_profile_model(params, exlist, reference)
      if predicted is None:
        predicted = self._predict_reflections(params, exlist)
        predicted = self._filter_reflections(params, exlist, predicted)
      if reference:
        predicted = self._match_with_reference(predicted, reference)
      extractor = self._create_extractor(params, exlist, predicted)

    # Initialise the integrator
    self._integrator = ReflectionBlockIntegrator(params, exlist, extractor)

  def integrate(self):
    ''' Integrate the reflections. '''
    return self._integrator.integrate()

  def _match_with_reference(self, predicted, reference):
    ''' Match predictions with reference spots. '''
    from dials.algorithms.peak_finding.spot_matcher import SpotMatcher
    from dials.util.command_line import Command
    Command.start("Matching reference spots with predicted reflections")
    match = SpotMatcher(max_separation=1)
    rind, pind = match(reference, predicted)
    h1 = predicted.select(pind)['miller_index']
    h2 = reference.select(rind)['miller_index']
    mask = rind == pind
    predicted.set_flags(pind.select(mask), reflections.flags.reference_spot)
    Command.end("Matched %d reference spots with predicted reflections" %
                mask.count(True))
    return predicted

  def _load_extractor(self, filename, params, exlist):
    ''' Load the shoebox extractor. '''
    from dials.model.serialize.reflection_block import ReflectionBlockExtractor
    assert(len(exlist) == 1)
    imageset = exlist[0].imageset
    return ReflectionBlockExtractor(
      filename,
      params.integration.shoebox.n_blocks,
      imageset)

  def _create_extractor(self, params, exlist, predicted):
    ''' Create the extractor. '''
    from dials.model.serialize.reflection_block import ReflectionBlockExtractor
    assert(len(exlist) == 1)
    imageset = exlist[0].imageset
    return ReflectionBlockExtractor(
      "shoebox.dat",
      params.integration.shoebox.n_blocks,
      imageset,
      predicted)

  def _compute_profile_model(self, params, experiments, reference):
    ''' Compute the profile model. '''
    from dials.algorithms.profile_model.profile_model import ProfileModel
    from math import pi
    if (params.integration.shoebox.sigma_b is None or
        params.integration.shoebox.sigma_m is None):
      assert(reference is not None)
      profile_model = ProfileModel(experiments[0], reference)
      params.integration.shoebox.sigma_b = profile_model.sigma_b() * 180.0 / pi
      params.integration.shoebox.sigma_m = profile_model.sigma_m() * 180.0 / pi
      print 'Sigma B: %f' % params.integration.shoebox.sigma_b
      print 'Sigma M: %f' % params.integration.shoebox.sigma_m

  def _predict_reflections(self, params, experiments):
    ''' Predict all the reflections. '''
    from dials.array_family import flex
    from math import pi
    n_sigma = params.integration.shoebox.n_sigma
    sigma_b = params.integration.shoebox.sigma_b * pi / 180.0
    sigma_m = params.integration.shoebox.sigma_m * pi / 180.0
    result = flex.reflection_table()
    for i, experiment in enumerate(experiments):
      predicted = flex.reflection_table.from_predictions(experiment)
      predicted['id'] = flex.size_t(len(predicted), i)
      predicted.compute_bbox(experiment, n_sigma, sigma_b, sigma_m)
      result.extend(predicted)
    return result

  def _filter_reflections(self, params, experiments, reflections):
    ''' Filter the reflections to integrate. '''
    from dials.util.command_line import Command
    from dials.algorithms import filtering
    from dials.array_family import flex

    # Set all reflections which overlap bad pixels to zero
    Command.start('Filtering reflections by detector mask')
    array_range = experiments[0].scan.get_array_range()
    mask = filtering.by_detector_mask(
      reflections['bbox'],
      experiments[0].imageset[0] >= 0,
      array_range)
    reflections.del_selected(mask != True)
    Command.end('Filtered %d reflections by detector mask' % len(reflections))

    # Filter the reflections by zeta
    min_zeta = params.integration.filter.by_zeta
    if min_zeta > 0:
      Command.start('Filtering reflections by zeta >= %f' % min_zeta)
      zeta = reflections.compute_zeta(experiments[0])
      reflections.del_selected(flex.abs(zeta) < min_zeta)
      n = len(reflections)
      Command.end('Filtered %d reflections by zeta >= %f' % (n, min_zeta))
      return reflections

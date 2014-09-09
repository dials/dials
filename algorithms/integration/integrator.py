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

  def __init__(self, params, experiments, profile_model, extractor=None):
    ''' Initialise the integrator. '''

    # Ensure we have 1 experiment at the moment
    assert(len(experiments) == 1)
    assert(extractor is not None)

    # Save the parameters
    self.params = params
    self.experiments = experiments
    self.extractor = extractor
    self.profile_model = profile_model

  def integrate(self):
    ''' Integrate all the reflections. '''
    from dials.array_family import flex
    from dials.algorithms.shoebox import MaskCode
    from dials.interfaces import BackgroundIface
    from dials.interfaces import IntensityIface
    from dials.interfaces import CentroidIface
    result = flex.reflection_table()
    flex.reflection_table._background_algorithm = flex.strategy(
      BackgroundIface.extension(self.params.integration.background.algorithm),
      self.params)
    flex.reflection_table._intensity_algorithm = flex.strategy(
      IntensityIface.extension(self.params.integration.intensity.algorithm),
      self.params)
    flex.reflection_table._centroid_algorithm = flex.strategy(
      CentroidIface.extension(self.params.integration.centroid.algorithm),
      self.params)
    for indices, reflections in self.extractor:
      reflections.compute_mask(self.experiments, self.profile_model)
      reflections.integrate(self.experiments, self.profile_model)
      bg_code = MaskCode.Valid | MaskCode.BackgroundUsed
      fg_code = MaskCode.Valid | MaskCode.Foreground
      n_bg = reflections['shoebox'].count_mask_values(bg_code)
      n_fg = reflections['shoebox'].count_mask_values(fg_code)
      reflections['n_background'] = n_bg
      reflections['n_foreground'] = n_fg
      del reflections['shoebox']
      del reflections['rs_shoebox']
      result.extend(reflections)
      print ''
    assert(len(result) > 0)
    result.sort('miller_index')
    result.compute_corrections(self.experiments)
    return result


class Integrator(object):
  ''' Integrate reflections '''

  def __init__(self, params, exlist, reference=None,
               predicted=None, shoeboxes=None):
    '''Initialise the script.'''
    from dials.algorithms.profile_model.profile_model import ProfileModelList

    # Load the reference spots and compute the profile parameters
    if reference:
      self._compute_profile_model(params, exlist, reference)
    else:
      self.profile_model = ProfileModelList.load(params)

    # Load the extractor based on the input
    if shoeboxes is not None:
      extractor = self._load_extractor(shoeboxes, params, exlist)
    else:
      if predicted is None:
        predicted = self._predict_reflections(params, exlist)
        predicted = self._filter_reflections(params, exlist, predicted)
      if reference:
        predicted = self._match_with_reference(predicted, reference)
      extractor = self._create_extractor(params, exlist, predicted)

    # Initialise the integrator
    self._integrator = ReflectionBlockIntegrator(
      params, exlist, self.profile_model, extractor)

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
    mask = (h1 == h2)
    predicted.set_flags(pind.select(mask), predicted.flags.reference_spot)
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
      params.integration.shoebox.block_size,
      imageset)

  def _create_extractor(self, params, exlist, predicted):
    ''' Create the extractor. '''
    from dials.model.serialize.reflection_block import ReflectionBlockExtractor
    assert(len(exlist) == 1)
    imageset = exlist[0].imageset
    return ReflectionBlockExtractor(
      "shoebox.dat",
      params.integration.shoebox.block_size,
      imageset,
      predicted)

  def _compute_profile_model(self, params, experiments, reference):
    ''' Compute the profile model. '''
    from dials.algorithms.profile_model.profile_model import ProfileModel
    self.profile_model = ProfileModel.compute(experiments[0], reference)
    for model in self.profile_model:
      print 'Sigma B: %f' % model.sigma_b(deg=True)
      print 'Sigma M: %f' % model.sigma_m(deg=True)

  def _predict_reflections(self, params, experiments):
    ''' Predict all the reflections. '''
    from dials.array_family import flex
    from dials.algorithms.profile_model.profile_model import ProfileModelList
    result = flex.reflection_table.from_predictions_multi(experiments)
    result.compute_bbox(experiments, self.profile_model)
    return result

  def _filter_reflections(self, params, experiments, reflections):
    ''' Filter the reflections to integrate. '''
    from dials.util.command_line import Command
    from dials.algorithms import filtering
    from dials.array_family import flex

    image = experiments[0].imageset[0]
    detector = experiments[0].detector
    if not isinstance(image, tuple):
      image = (image,)
    image_mask = []
    for im, panel in zip(image, detector):
      tr = panel.get_trusted_range()
      m = im > int(tr[0])
      image_mask.append(m)
    image_mask = tuple(image_mask)

    # Set all reflections which overlap bad pixels to zero
    Command.start('Filtering reflections by detector mask')
    array_range = experiments[0].scan.get_array_range()
    mask = filtering.by_detector_mask(
      reflections['panel'],
      reflections['bbox'],
      image_mask,
      array_range)
    reflections.del_selected(mask != True)
    Command.end('Filtered %d reflections by detector mask' % len(reflections))
    assert(len(reflections) > 0)

    # Filter the reflections by zeta
    min_zeta = params.integration.filter.by_zeta
    if min_zeta > 0:
      Command.start('Filtering reflections by zeta >= %f' % min_zeta)
      zeta = reflections.compute_zeta(experiments[0])
      reflections.del_selected(flex.abs(zeta) < min_zeta)
      n = len(reflections)
      Command.end('Filtered %d reflections by zeta >= %f' % (n, min_zeta))
      assert(len(reflections) > 0)
    return reflections

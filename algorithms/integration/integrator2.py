

from __future__ import division

class Integrator(object):
  ''' A class to perform the integration. '''

  def __init__(self, params, experiments, reference=None, extractor=None):
    ''' Initialise the integrator. '''

    # Ensure we have 1 experiment at the moment
    assert(len(experiments) == 1)

    # Save the parameters
    self.params = params

    # Filter the reference profiles
    reference = self.filter_reference(reference)

    # If no extractor is given, then compute the profile model and predict the
    # reflections and filter them.
    if extractor is None:

      # Compute the profile model
      self.compute_profile_model(params, experiments, reference)

      # Predict the reflections and compute the bounding box
      predictions = self.predict_reflections(params, experiments)

      # Filter the reflections
      predictions = self.filter_reflections(params, experiments, predictions)

      # Initialise the extractor
      self.extractor = []

  def integrate(self):
    ''' Integrate all the reflections. '''
    from dials.array_family import flex
    result = flex.reflection_table()
    for reflections in self.extractor:
      reflections.integrate(self.experiments[0])
      del reflections['shoebox']
      del reflections['rs_shoebox']
      result.extend(reflections)
    result.sort('miller_index')
    return result

  def filter_reference(self, reference):
    ''' Filter the reference profiles. '''
    from dials.array_family import flex
    from dials.util.command_line import Command
    Command.start('Removing invalid coordinates')
    xyz = reference['xyzcal.mm']
    mask = flex.bool([x == (0, 0, 0) for x in xyz])
    reference.del_selected(mask)
    Command.end('Removed invalid coordinates, %d remaining' % len(reference))
    return reference

  def compute_profile_model(self, params, experiments, reference):
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

  def predict_reflections(self, params, experiments):
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

  def filter_reflections(self, params, experiments, reflections):
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

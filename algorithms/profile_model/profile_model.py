#
# profile_model.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division
from libtbx.phil import parse

phil_scope = parse('''
  profile
    .multiple=True
    .optional=False
  {
    n_sigma = 3
      .help = "The number of standard deviations of the beam divergence and the"
              "mosaicity to use for the bounding box size."
      .type = float

    sigma_b = 0
      .help = "The E.S.D. of the beam divergence"
      .type = float

    sigma_m = 0
      .help = "The E.S.D. of the reflecting range"
      .type = float
  }
''')

class ProfileModel(object):
  ''' A class to encapsulate the profile model. '''

  def __init__(self, n_sigma, sigma_b, sigma_m):
    ''' Initialise with the parameters. '''
    self._n_sigma = n_sigma
    self._sigma_b = sigma_b
    self._sigma_m = sigma_m
    assert(self._n_sigma > 0)
    assert(self._sigma_b > 0)
    assert(self._sigma_m > 0)

  def sigma_b(self, deg=True):
    ''' Return sigma_b. '''
    from math import pi
    if deg == True:
      return self._sigma_b * 180.0 / pi
    return self._sigma_b

  def sigma_m(self, deg=True):
    ''' Return sigma_m. '''
    from math import pi
    if deg == True:
      return self._sigma_m * 180.0 / pi
    return self._sigma_m

  def n_sigma(self):
    ''' The number of sigmas. '''
    return self._n_sigma

  def compute_bbox(self, experiment, reflections, sigma_b_multiplier=2.0):
    ''' Compute the bounding box. '''
    from dials.algorithms.shoebox import BBoxCalculator

    # Check the input
    assert(sigma_b_multiplier >= 1.0)

    # Compute the size in reciprocal space. Add a sigma_b multiplier to enlarge
    # the region of background in the shoebox
    delta_b = self._n_sigma * self._sigma_b * sigma_b_multiplier
    delta_m = self._n_sigma * self._sigma_m

    # Create the bbox calculator
    calculate = BBoxCalculator(
      experiment.beam,
      experiment.detector,
      experiment.goniometer,
      experiment.scan,
      delta_b,
      delta_m)

    # Calculate the bounding boxes of all the reflections
    reflections['bbox'] = calculate(
      reflections['s1'],
      reflections['xyzcal.mm'].parts()[2],
      reflections['panel'])

  def compute_partiality(self, experiment, reflections):
    ''' Compute the partiality. '''
    from dials.algorithms.shoebox import PartialityCalculator

    # Compute the size in reciprocal space.
    delta_m = self._n_sigma * self._sigma_m

    # Create the partiality calculator
    calculate = PartialityCalculator(
      experiment.beam,
      experiment.goniometer,
      experiment.scan,
      delta_m)

    # Compute the partiality
    reflections['partiality'] = calculate(
      reflections['s1'],
      reflections['xyzcal.mm'].parts()[2],
      reflections['bbox'])

  def compute_mask(self, experiment, reflections):
    ''' Compute the shoebox mask. '''
    from dials.algorithms.shoebox import MaskForeground

    # Compute the size in reciprocal space. Add a sigma_b multiplier to enlarge
    # the region of background in the shoebox
    delta_b = self._n_sigma * self._sigma_b
    delta_m = self._n_sigma * self._sigma_m

    # Create the mask calculator
    mask_foreground = MaskForeground(
      experiment.beam,
      experiment.detector,
      experiment.goniometer,
      experiment.scan,
      delta_d,
      delta_m)

    # Mask the foreground region
    mask_foreground(
      reflections['shoebox'],
      reflections['s1'],
      reflections['xyzcal.px'].parts()[2],
      reflections['panel'])

  @classmethod
  def compute(cls, experiment, reflections, min_zeta=0.05):
    ''' Compute the profile model. '''
    from dials.algorithms.profile_model.profile_model_calculator \
      import ProfileModelCalculator
    calculator = ProfileModelCalculator(experiment, reflections, min_zeta)
    n_sigma = 3
    sigma_b = calculator.sigma_b()
    sigma_m = calculator.sigma_m()
    return cls(n_sigma, sigma_b, sigma_m)


class ProfileModelList(object):
  ''' A class to represent multiple profile models. '''

  def __init__(self):
    ''' Initialise the model list. '''
    self._models = []

  def __getitem__(self, index):
    ''' Get a profile model. '''
    return self._models[index]

  def __len__(self):
    ''' Return the number of models. '''
    return len(self._models)

  def __iter__(self):
    ''' Iterate through the models. '''
    for i in range(len(self)):
      yield self[i]

  def append(self, model):
    ''' Add another model. '''
    self._models.append(model)

  def compute_bbox(self, experiments, reflections, sigma_b_multiplier=2.0):
    ''' Compute the bounding boxes. '''
    from dials.algorithms.shoebox import BBoxMultiCalculator
    from dials.algorithms.shoebox import BBoxCalculator

    # Check the input
    assert(sigma_b_multiplier >= 1.0)
    assert(len(experiments) == len(self))

    # The class to do the calculation
    calculate = BBoxMultiCalculator()

    # Loop through the experiments and models
    for experiment, model in zip(experiments, self):

      # Compute the size in reciprocal space. Add a sigma_b multiplier to enlarge
      # the region of background in the shoebox
      delta_b = model.n_sigma() * model.sigma_b(deg=False) * sigma_b_multiplier
      delta_m = model.n_sigma() * model.sigma_m(deg=False)

      # Create the bbox calculator
      calculate.append(BBoxCalculator(
        experiment.beam,
        experiment.detector,
        experiment.goniometer,
        experiment.scan,
        delta_b,
        delta_m))

    # Calculate the bounding boxes of all the reflections
    reflections['bbox'] = calculate(
      reflections['id'],
      reflections['s1'],
      reflections['xyzcal.mm'].parts()[2],
      reflections['panel'])

  def compute_partiality(self, experiments, reflections):
    ''' Compute the partiality. '''
    from dials.algorithms.shoebox import PartialityMultiCalculator
    from dials.algorithms.shoebox import PartialityCalculator

    # Check the input
    assert(len(experiments) == len(self))

    # The partiality calculator
    calculate = PartialityMultiCalculator()

    # Loop through the experiments and models
    for experiment, model in zip(experiments, self):

      # Compute the size in reciprocal space.
      delta_m = model.n_sigma() * model.sigma_m(deg=False)

      # Create the partiality calculator
      calculate.append(PartialityCalculator(
        experiment.beam,
        experiment.goniometer,
        experiment.scan,
        delta_m))

    # Compute the partiality
    reflections['partiality'] = calculate(
      reflections['id'],
      reflections['s1'],
      reflections['xyzcal.mm'].parts()[2],
      reflections['bbox'])

  def compute_mask(self, experiments, reflections):
    ''' Compute the shoebox mask. '''
    from dials.algorithms.shoebox import MaskMultiForeground
    from dials.algorithms.shoebox import MaskForeground

    # Check the input
    assert(len(experiments) == len(self))

    # The class to do the calculation
    mask_foreground = MaskMultiForeground()

    # Loop through the experiments and models
    for experiment, model in zip(experiments, self):

      # Compute the size in reciprocal space.
      delta_b = model.n_sigma() * model.sigma_b(deg=False)
      delta_m = model.n_sigma() * model.sigma_m(deg=False)

      # Create the mask calculator
      mask_foreground.append(MaskForeground(
        experiment.beam,
        experiment.detector,
        experiment.goniometer,
        experiment.scan,
        delta_b,
        delta_m))

    # Mask the foreground region
    mask_foreground(
      reflections['id'],
      reflections['shoebox'],
      reflections['s1'],
      reflections['xyzcal.px'].parts()[2],
      reflections['panel'])

  @classmethod
  def compute(cls, experiments, reflections, min_zeta=0.05):
    ''' Compute the profile models. '''
    from dials.util.command_line import heading
    assert(len(experiments) > 0)

    print "=" * 80
    print ""
    print heading("Computing Profile Model")
    print ""

    # Split the reflections by experiment id
    if len(experiments) > 1:
      reflections_split = reflections.split_by_experiment_id()
      assert(len(reflections_split) == len(experiments))
    else:
      reflections_split = [reflections]

    # Compute the profile models
    profile_models = cls()
    for exp, ref in zip(experiments, reflections_split):
      profile_models.append(ProfileModel.compute(exp, ref, min_zeta))

    # Return the profile models
    return profile_models

  @classmethod
  def load(cls, params):
    ''' Load from phil parameters. '''
    from math import pi
    assert(len(params.profile) > 1)
    profile_model = cls()
    for i in range(1, len(params.profile)):
      profile_model.append(ProfileModel(
        params.profile[i].n_sigma,
        params.profile[i].sigma_b * pi / 180.0,
        params.profile[i].sigma_m * pi / 180.0))
    return profile_model

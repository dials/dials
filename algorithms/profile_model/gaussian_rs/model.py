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
from dials.algorithms.profile_model.interface import ProfileModelIface

phil_scope = parse('''

  gaussian_rs
  {
    filter
    {
      min_zeta = 0.05
        .type = float
        .help = "Filter reflections by min zeta"
    }

    model
      .multiple = True
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
  }

''')

class ProfileModel(object):
  ''' A class to encapsulate the profile model. '''

  def __init__(self, n_sigma, sigma_b, sigma_m, deg=False):
    ''' Initialise with the parameters. '''
    from math import pi
    self._n_sigma = n_sigma
    if deg == True:
      self._sigma_b = sigma_b * pi / 180.0
      self._sigma_m = sigma_m * pi / 180.0
    else:
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

  def delta_b(self, deg=True):
    ''' Return delta_b. '''
    return self.sigma_b(deg) * self.n_sigma()

  def delta_m(self, deg=True):
    ''' Return delta_m. '''
    return self.sigma_m(deg) * self.n_sigma()

  def compute_bbox(self, experiment, reflections, sigma_b_multiplier=2.0):
    ''' Compute the bounding box. '''
    from dials.algorithms.profile_model.gaussian_rs import BBoxCalculator

    # Check the input
    assert(sigma_b_multiplier >= 1.0)

    # Compute the size in reciprocal space. Add a sigma_b multiplier to enlarge
    # the region of background in the shoebox
    delta_b = self._n_sigma * self._sigma_b * sigma_b_multiplier
    delta_m = self._n_sigma * self._sigma_m

    # Create the bbox calculator
    calculate = BBoxCalculator(experiment, delta_b, delta_m)

    # Calculate the bounding boxes of all the reflections
    reflections['bbox'] = calculate(
      reflections['s1'],
      reflections['xyzcal.px'].parts()[2],
      reflections['panel'])

  def compute_partiality(self, experiment, reflections):
    ''' Compute the partiality. '''
    from dials.algorithms.profile_model.gaussian_rs import PartialityCalculator

    # Create the partiality calculator
    calculate = PartialityCalculator(experiment, sigma_m(deg=False))

    # Compute the partiality
    reflections['partiality'] = calculate(
      reflections['s1'],
      reflections['xyzcal.px'].parts()[2],
      reflections['bbox'])

  def compute_mask(self, experiment, reflections):
    ''' Compute the shoebox mask. '''
    from dials.algorithms.profile_model.gaussian_rs import MaskCalculator

    # Compute the size in reciprocal space. Add a sigma_b multiplier to enlarge
    # the region of background in the shoebox
    delta_b = self._n_sigma * self._sigma_b
    delta_m = self._n_sigma * self._sigma_m

    # Create the mask calculator
    mask_foreground = MaskCalculator(experiment, delta_b, delta_m)

    # Mask the foreground region
    mask_foreground(
      reflections['shoebox'],
      reflections['s1'],
      reflections['xyzcal.px'].parts()[2],
      reflections['panel'])

  @classmethod
  def compute(cls, experiment, reflections, min_zeta=0.05):
    ''' Compute the profile model. '''
    from dials.algorithms.profile_model.gaussian_rs.calculator \
      import ProfileModelCalculator
    calculator = ProfileModelCalculator(experiment, reflections, min_zeta)
    n_sigma = 3
    sigma_b = calculator.sigma_b()
    sigma_m = calculator.sigma_m()
    return cls(n_sigma, sigma_b, sigma_m)


class ProfileModelList(ProfileModelIface):
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

  def predict_reflections(self, experiments, dmin=None, dmax=None, margin=1,
                          force_static=False):
    ''' Predict the reflections. '''
    from dials.array_family import flex
    return flex.reflection_table.from_predictions_multi(
      experiments,
      dmin=dmin,
      dmax=dmax,
      margin=margin,
      force_static=force_static)

  def compute_bbox(self, experiments, reflections, sigma_b_multiplier=2.0):
    ''' Compute the bounding boxes. '''
    from dials.algorithms.profile_model.gaussian_rs import BBoxMultiCalculator
    from dials.algorithms.profile_model.gaussian_rs import BBoxCalculator

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
      calculate.append(BBoxCalculator(experiment, delta_b, delta_m))

    # Calculate the bounding boxes of all the reflections
    reflections['bbox'] = calculate(
      reflections['id'],
      reflections['s1'],
      reflections['xyzcal.px'].parts()[2],
      reflections['panel'])

  def compute_partiality(self, experiments, reflections):
    ''' Compute the partiality. '''
    from dials.algorithms.profile_model.gaussian_rs import PartialityMultiCalculator
    from dials.algorithms.profile_model.gaussian_rs import PartialityCalculator

    # Check the input
    assert(len(experiments) == len(self))

    # The partiality calculator
    calculate = PartialityMultiCalculator()

    # Loop through the experiments and models
    # Create the partiality calculator
    for experiment, model in zip(experiments, self):
      calculate.append(PartialityCalculator(experiment, model.sigma_m(deg=False)))

    # Compute the partiality
    reflections['partiality'] = calculate(
      reflections['id'],
      reflections['s1'],
      reflections['xyzcal.px'].parts()[2],
      reflections['bbox'])

  def compute_mask(self, experiments, reflections):
    ''' Compute the shoebox mask. '''
    from dials.algorithms.profile_model.gaussian_rs import MaskMultiCalculator
    from dials.algorithms.profile_model.gaussian_rs import MaskCalculator

    # Check the input
    assert(len(experiments) == len(self))

    # The class to do the calculation
    mask_foreground = MaskMultiCalculator()

    # Loop through the experiments and models
    for experiment, model in zip(experiments, self):

      # Compute the size in reciprocal space.
      delta_b = model.n_sigma() * model.sigma_b(deg=False)
      delta_m = model.n_sigma() * model.sigma_m(deg=False)

      # Create the mask calculator
      mask_foreground.append(MaskCalculator(experiment, delta_b, delta_m))

    # Mask the foreground region
    mask_foreground(
      reflections['id'],
      reflections['shoebox'],
      reflections['s1'],
      reflections['xyzcal.px'].parts()[2],
      reflections['panel'])

  @classmethod
  def compute(cls, params, experiments, reflections):
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
    min_zeta = params.gaussian_rs.filter.min_zeta
    profile_models = cls()
    for exp, ref in zip(experiments, reflections_split):
      model = ProfileModel.compute(exp, ref, min_zeta)
      print " Sigma_b: %.3f degrees" % model.sigma_b(deg=True)
      print " Sigma_m: %.3f degrees" % model.sigma_m(deg=True)
      profile_models.append(model)

    # Return the profile models
    return profile_models

  @classmethod
  def load(cls, params):
    ''' Load from phil parameters. '''
    from math import pi
    assert(len(params.gaussian_rs.model) > 0)
    profile_model = cls()
    for i in range(len(params.gaussian_rs.model)):
      profile_model.append(ProfileModel(
        params.gaussian_rs.model[i].n_sigma,
        params.gaussian_rs.model[i].sigma_b * pi / 180.0,
        params.gaussian_rs.model[i].sigma_m * pi / 180.0))
    return profile_model

  @classmethod
  def create(cls, params, experiments, reflections=None):
    ''' Create the profile model. '''
    if len(params.profile.gaussian_rs.model) > 0:
      assert(len(params.profile.gaussian_rs.model) == len(experiments))
      model = ProfileModelList.load(params.profile)
    else:
      assert(reflections is not None)
      model = ProfileModelList.compute(params.profile, experiments, reflections)
    return model

  def dump(self):
    ''' Dump the profile model to phil parameters. '''
    from dials.algorithms.profile_model import factory
    phil_str = '\n'.join([
      '''
      profile {
        gaussian_rs {
          model {
            n_sigma=%g
            sigma_b=%g
            sigma_m=%g
          }
        }
      }
      ''' % (
        m.n_sigma(),
        m.sigma_b(deg=True),
        m.sigma_m(deg=True)) for m in self])
    return factory.phil_scope.fetch(source=parse(phil_str))

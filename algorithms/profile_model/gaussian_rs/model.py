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
    scan_varying = False
        .type = bool
        .help = "Calculate a scan varying model"

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

    scan_varying_model
      .multiple = True
    {
      n_sigma = 3
        .help = "The number of standard deviations of the beam divergence and the"
                "mosaicity to use for the bounding box size."
        .type = float

      sigma_b = None
        .help = "The E.S.D. of the beam divergence"
        .type = floats

      sigma_m = None
        .help = "The E.S.D. of the reflecting range"
        .type = floats
    }
  }

''')

class ProfileModel(ProfileModelIface):
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

  def predict_reflections(self, experiment, dmin=None, dmax=None, margin=1,
                          force_static=False, **kwargs):
    ''' Predict the reflections. '''
    from dials.array_family import flex
    return flex.reflection_table.from_predictions(
      experiment,
      dmin=dmin,
      dmax=dmax,
      margin=margin,
      force_static=force_static)

  def compute_bbox(self, experiment, reflections,
                   sigma_b_multiplier=2.0, **kwargs):
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
    bbox = calculate(
      reflections['s1'],
      reflections['xyzcal.px'].parts()[2],
      reflections['panel'])

    # Return the bounding boxes
    return bbox

  def compute_partiality(self, experiment, reflections, **kwargs):
    ''' Compute the partiality. '''
    from dials.algorithms.profile_model.gaussian_rs import PartialityCalculator

    # Create the partiality calculator
    calculate = PartialityCalculator(experiment, self.sigma_m(deg=False))

    # Compute the partiality
    partiality = calculate(
      reflections['s1'],
      reflections['xyzcal.px'].parts()[2],
      reflections['bbox'])

    # Return the partiality
    return partiality

  def compute_mask(self, experiment, reflections, **kwargs):
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

  def dump(self):
    ''' Dump the profile model to phil parameters. '''
    from dials.algorithms.profile_model import factory
    phil_str = '''
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
        self.n_sigma(),
        self.sigma_b(deg=True),
        self.sigma_m(deg=True))
    return factory.phil_scope.fetch(source=parse(phil_str))


class ScanVaryingProfileModel(ProfileModelIface):
  ''' A class to encapsulate the profile model. '''

  def __init__(self, n_sigma, sigma_b, sigma_m, deg=False, num_used=None):
    ''' Initialise with the parameters. '''
    from math import pi
    self._num_used = num_used
    self._n_sigma = n_sigma
    if deg == True:
      self._sigma_b = sigma_b * pi / 180.0
      self._sigma_m = sigma_m * pi / 180.0
    else:
      self._sigma_b = sigma_b
      self._sigma_m = sigma_m
    assert(self._n_sigma > 0)
    assert(self._sigma_b.all_gt(0))
    assert(self._sigma_m.all_gt(0))
    assert(len(self._sigma_b) == len(self._sigma_m))

  def sigma_b(self, index=None, deg=True):
    ''' Return sigma_b. '''
    from math import pi
    if index is None:
      sigma_b = self._sigma_b
    else:
      sigma_b = self._sigma_b[index]
    if deg == True:
      return sigma_b * 180.0 / pi
    return sigma_b

  def sigma_m(self, index=None, deg=True):
    ''' Return sigma_m. '''
    from math import pi
    if index is None:
      sigma_m = self._sigma_m
    else:
      sigma_m = self._sigma_m[index]
    if deg == True:
      return sigma_m * 180.0 / pi
    return sigma_m

  def n_sigma(self):
    ''' The number of sigmas. '''
    return self._n_sigma

  def delta_b(self, index=None, deg=True):
    ''' Return delta_b. '''
    return self.sigma_b(index, deg=deg) * self.n_sigma()

  def delta_m(self, index=None, deg=True):
    ''' Return delta_m. '''
    return self.sigma_m(index, deg=deg) * self.n_sigma()

  def num_used(self, index=None):
    ''' Return number of reflections used. '''
    if index is None:
      return self._num_used
    return self._num_used[index]

  def __len__(self):
    assert(len(self._sigma_m) == len(self._sigma_b))
    return len(self._sigma_m)

  def predict_reflections(self, experiment, dmin=None, dmax=None, margin=1,
                          force_static=False, **kwargs):
    ''' Predict the reflections. '''
    from dials.array_family import flex
    return flex.reflection_table.from_predictions(
      experiment,
      dmin=dmin,
      dmax=dmax,
      margin=margin,
      force_static=force_static)

  def compute_bbox(self, experiment, reflections,
                   sigma_b_multiplier=2.0, **kwargs):
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
    bbox = calculate(
      reflections['s1'],
      reflections['xyzcal.px'].parts()[2],
      reflections['panel'])

    # Return the bounding boxes
    return bbox

  def compute_partiality(self, experiment, reflections, **kwargs):
    ''' Compute the partiality. '''
    from dials.algorithms.profile_model.gaussian_rs import PartialityCalculator

    # Create the partiality calculator
    calculate = PartialityCalculator(experiment, self._sigma_m)

    # Compute the partiality
    partiality = calculate(
      reflections['s1'],
      reflections['xyzcal.px'].parts()[2],
      reflections['bbox'])

    # Return the partiality
    return partiality

  def compute_mask(self, experiment, reflections, **kwargs):
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

  def dump(self):
    ''' Dump the profile model to phil parameters. '''
    from dials.algorithms.profile_model import factory
    phil_str = '''
      profile {
        gaussian_rs {
          scan_varying=True
          scan_varying_model {
            n_sigma=%s
            sigma_b=%s
            sigma_m=%s
          }
        }
      }
      ''' % (
        self.n_sigma(),
        ','.join(["%g" % v for v in self.sigma_b(deg=True)]),
        ','.join(["%g" % v for v in self.sigma_m(deg=True)]))
    return factory.phil_scope.fetch(source=parse(phil_str))

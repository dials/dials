#!/usr/bin/env python
#
# summation_integration_ext.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import division

from dials.interfaces import IntensityIface

class SummationIntegrationExt(IntensityIface):
  ''' Extension class to provide 3d summation integration.

  Summation integration should be performed anyway so this extension is just a
  placeholder to give some options for when we only want to do summation.

  '''

  name = 'sum'

  @classmethod
  def phil(cls):
    from libtbx.phil import parse
    phil = parse('''

      integrator = *auto 3d flat3d 2d single2d stills
        .type = choice
        .help = "The integrator to use."
        .expert_level=3

    ''')
    return phil

  def __init__(self, params, experiments, profile_model):
    '''
    Initialise the algorithm.

    :param params: The input parameters
    :param experiments: The list of experiments
    :param profile_model: The profile model

    '''
    pass

  def compute_intensity(self, reflections):
    '''
    Compute the intensity.

    :param reflections: The list of reflections

    '''
    pass

  @classmethod
  def type(cls, params, experiments):
    '''
    Return the type of the integrator.

    :param params: The input parameters
    :param experiments: The list of experiments

    '''
    from libtbx import Auto
    integrator_type = params.integration.intensity.sum.integrator
    if integrator_type == Auto or integrator_type == 'auto':
      if experiments.all_stills():
        integrator_type = 'stills'
      elif experiments.all_sweeps():
        integrator_type = '3d'
      else:
        raise RuntimeError("Experiments must be all sweeps or all stills.")
    return integrator_type

#!/usr/bin/env cctbx.python
#
#  Copyright (C) (2016) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman.
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#

"""Contains classes used to manage the reflections used during scaling,
principally ObservationManager."""

from __future__ import division
from dials.array_family import flex
import cctbx.crystal
import cctbx.miller
from dials_scaling_helpers_ext import (GroupedObservations,
  minimum_multiplicity_selection)
from libtbx import phil
from libtbx.utils import Sorry

observation_manager_phil_str='''
observations
{
  integration_type = *sum prf mix
    .type = choice

  i_over_sigma_cutoff = 2
    .type = float(value_min=0)

  min_multiplicity = 3
    .type = int(value_min=0)
}
'''

om_scope = phil.parse(observation_manager_phil_str)

class ObservationManager(object):
  """A class to maintain information about the reflections during
  optimisation of their scale factors."""

  def __init__(self, reflections, experiment, params=None):

    if params is None: params = om_scope.extract()

    # initial filter
    reflections = reflections.select(reflections.get_flags(
      reflections.flags.integrated))

    # create new column containing the reduced Miller index
    xl = experiment.crystal
    symm = cctbx.crystal.symmetry(xl.get_unit_cell(),
                                  space_group=xl.get_space_group())
    hkl_set = cctbx.miller.set(symm, reflections['miller_index'])
    asu_set = hkl_set.map_to_asu()
    reflections['asu_miller_index'] = asu_set.indices()

    # sort table by reduced hkl
    reflections.sort('asu_miller_index')

    if params.observations.integration_type == 'mix':
      raise Sorry('integration_type=mix is not supported yet')
    else:
      ikey = 'intensity.' + params.observations.integration_type + '.value'
      vkey = 'intensity.' + params.observations.integration_type + '.variance'

    # filters
    sel = reflections[vkey] > 0
    reflections = reflections.select(sel)
    if params.observations.i_over_sigma_cutoff > 0:
      ios = reflections[ikey] / flex.sqrt(reflections[vkey])
      sel = ios >= params.observations.i_over_sigma_cutoff
      reflections = reflections.select(sel)
    if params.observations.min_multiplicity > 0:
      sel = minimum_multiplicity_selection(reflections['asu_miller_index'],
        params.observations.min_multiplicity)
      reflections = reflections.select(sel)

    # extract columns of interest
    gp_idx = reflections['asu_miller_index']
    intensity = reflections[ikey]
    weight = 1. / reflections[vkey]
    phi = reflections['xyzcal.mm'].parts()[2]
    scale = flex.double(len(reflections), 1.0)

    # set up reflection grouping object
    self._go = GroupedObservations(gp_idx, intensity, weight, phi, scale)

    return

  @property
  def intensity(self):
    return self._go.get_intensity()

  @intensity.setter
  def intensity(self, intensity):
    self._go.set_intensity(intensity)

  @property
  def weight(self):
    return self._go.get_weight()

  @weight.setter
  def weight(self, weight):
    self._go.set_weight(weight)

  @property
  def scale(self):
    return self._go.get_scale()

  @scale.setter
  def scale(self, scale):
    self._go.set_scale(scale)

  @property
  def group_index(self):
    return self._go.get_group_index()

  @property
  def phi(self):
    return self._go.get_phi()

  @property
  def group_size(self):
    return self._go.get_group_size()

  def get_average_intensity(self):
    '''Calculate the weighted average intensity in reflection groups at the
    current scales using the formula of Hamilton, Rollett and Sparks (1965).
    Return as a vector for each observation'''

    # delegate to the GroupedObservations object
    return self._go.get_average_intensity()


class Script(object):
  '''A class for running the script.'''

  def __init__(self):
    '''Initialise the script.'''
    from dials.util.options import OptionParser
    from libtbx.phil import parse
    import libtbx.load_env

    # The script usage
    import __main__
    usage  = ("usage: dials.python {0} [options] [param.phil] "
              "integrated.pickle refined_experiments.json").format(
      __main__.__file__)

    # Create the parser
    self.parser = OptionParser(
      usage=usage,
      #phil=working_phil,
      read_reflections=True,
      read_experiments=True,
      check_format=False,
      epilog="Test script for the ObservationManager class.")

    return

  def run(self):
    '''Execute the script.'''

    from dials.util.options import flatten_reflections, flatten_experiments

    # Parse the command line
    params, options = self.parser.parse_args(show_diff_phil=False)
    reflections = flatten_reflections(params.input.reflections)
    experiments = flatten_experiments(params.input.experiments)

    # Try to load the models and data
    nexp = len(experiments)
    if nexp == 0:
      print "No Experiments found in the input"
      self.parser.print_help()
      return
    if len(reflections) == 0:
      print "No reflection data found in the input"
      self.parser.print_help()
      return
    if len(reflections) > 1:
      raise Sorry("Only one reflections list can be imported at present")
    reflections = reflections[0] # first reflection list
    reflections = reflections.select(reflections['id'] == 0) # first experiment
    if len(reflections) == 0:
      print "No reflection data for the first experiment found in the input"
      self.parser.print_help()
      return
    if len(experiments) > 1:
      raise Sorry("Only one experiment can be imported at present")
    experiment = experiments[0]

    om = ObservationManager(reflections, experiment)

    om.group_index

# For testing, instantiate from reflections passed at the command line.
if __name__ == '__main__':
  from dials.util import halraiser
  try:
    script = Script()
    script.run()
  except Exception as e:
    halraiser(e)

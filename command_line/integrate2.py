#!/usr/bin/env python
#
# dials.integrate.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division
from dials.util.script import ScriptRunner

class Script(ScriptRunner):
  ''' The integration program. '''

  def __init__(self):
    '''Initialise the script.'''

    # The script usage
    usage  = "usage: %prog [options] experiment.json"

    # Initialise the base class
    ScriptRunner.__init__(self, usage=usage)

    # Output filename option
    self.config().add_option(
      '-o', '--output',
      dest = 'output',
      type = 'string', default = 'integrated.pickle',
      help = 'Set the filename for integrated reflections.')

    # The predicted reflections to integrate
    self.config().add_option(
      '-p', '--predicted',
      dest = 'predicted',
      type = 'string', default = None,
      help = 'Specify predicted reflections.')

    # The intermediate shoebox data
    self.config().add_option(
      '-r', '--reference',
      dest = 'reference',
      type = 'string', default = None,
      help = 'Specify reference reflections.')

    # The intermediate shoebox data
    self.config().add_option(
      '-s', '--shoeboxes',
      dest = 'shoeboxes',
      type = 'string', default = None,
      help = 'Specify shoeboxes to integrate.')

  def main(self, params, options, args):
    ''' Perform the integration. '''
    from dials.algorithms.integration.integrator2 import Integrator
    from time import time

    # Check the number of arguments is correct
    start_time = time()

    # Check the number of command line arguments
    if len(args) != 1:
      self.config().print_help()

    # Load the experiment list
    exlist = self.load_experiments(args[0])

    # Load the extractor
    extractor = self.extractor(params, options, exlist)

    # Initialise the integrator
    integrator = Integrator(params, exlist, extractor)

    # Integrate the reflections
    reflections = integrator.integrate()

    # Save the reflections
    self.save_reflections(reflections, options.output)

    # Print the total time taken
    print "\nTotal time taken: ", time() - start_time

  def extractor(self, params, options, exlist):
    ''' Get the extractor. '''
    if options.shoeboxes:
      extractor = self.load_extractor(options.shoeboxes)
    else:
      reference = self.load_reference(options.reference)
      if reference:
        self.compute_profile_model(params, exlist, reference)
      if options.predicted:
        predicted = self.load_predicted(options.predicted)
      else:
        predicted = self.predict_reflections(params, exlist)
        predicted = self.filter_reflections(params, exlist, predicted)
      extractor = self.create_extractor(exlist, predicted)
    return extractor

  def load_predicted(self, filename):
    ''' Load the predicted reflections. '''
    from dials.array_family import flex
    return flex.reflection_table.from_pickle(filename)

  def load_experiments(self, filename):
    ''' Load the experiment list. '''
    from dials.model.experiment.experiment_list import ExperimentListFactory
    from dials.util.command_line import Command
    Command.start('Loading experiments from %s' % filename)
    exlist = ExperimentListFactory.from_json_file(filename)
    Command.end('Loaded experiments from %s' % filename)
    if len(exlist) == 0:
      raise RuntimeError('experiment list is empty')
    elif len(exlist.imagesets()) > 1 or len(exlist.detectors()) > 1:
      raise RuntimeError('experiment list contains > 1 imageset or detector')
    return exlist

  def load_reference(self, filename):
    ''' Load the reference spots. '''
    from dials.util.command_line import Command
    from dials.array_family import flex
    if filename:
      Command.start('Loading reference spots from %s' % filename)
      reference = flex.reflection_table.from_pickle(filename)
      Command.end('Loaded reference spots from %s' % filename)
      Command.start('Removing reference spots with invalid coordinates')
      xyz = reference['xyzcal.mm']
      mask = flex.bool([x == (0, 0, 0) for x in xyz])
      reference.del_selected(mask)
      Command.end('Removed reference spots with invalid coordinates, \
                  %d remaining' % len(reference))
    else:
      reference = None
    return reference

  def load_extractor(self, filename):
    ''' Load the shoebox extractor. '''
    return None

  def create_extractor(self, exlist, predicted):
    ''' Create the extractor. '''
    from dials.model.serialize.reflection_block import ReflectionBlockExtractor
    assert(len(exlist) == 1)
    imageset = exlist[0].imageset
    return ReflectionBlockExtractor(
      "shoebox.dat",
      imageset,
      predicted,
      1)

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

  def save_reflections(self, reflections, filename):
    ''' Save the reflections to file. '''
    from dials.util.command_line import Command
    Command.start('Saving %d reflections to %s' % (len(reflections), filename))
    reflections.as_pickle(filename)
    Command.end('Saved %d reflections to %s' % (len(reflections), filename))


if __name__ == '__main__':
  script = Script()
  script.run()

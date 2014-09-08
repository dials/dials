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

help_message = '''

This program is used to integrate the reflections on the diffraction images. It
is called with an experiment list outputted from dials.index or dials.refine.
The extend of the shoeboxes is specified through the profile parameters
shoebox.sigma_b and shoebox.sigma_m (use the --show-config option for more
details). These parameters can be specified directly, otherwise a set of strong
indexed reflections are needed to form the profile model; these are specified
using the -r (for reference) option. The program can also be called with a
specific set of predictions using the -p option.

Once a profile model is given and the size of the measurement boxes have been
calculated, the program will extract the reflections to file. The reflections
will then be integrated. The reflections can be integrated with different
options using the same measurement boxes by giving the measurement box file
using the -s option. This will skip reading the measurement boxes and go
directly to integrating the reflections.

Examples:

  dials.integrate experiments.json -r indexed.pickle

  dials.integrate experiments.json -r indexed.pickle -o integrated.pickle

  dials.integrate experiments.json shoebox.sigma_b=0.024 shoebox.sigma_m=0.044

  dials.integrate experiments.json -p predicted.pickle -r indexed.pickle

  dials.integrate experiments.json -s shoeboxes.dat

'''


class Script(object):
  ''' The integration program. '''

  def __init__(self):
    '''Initialise the script.'''
    from dials.util.options import OptionParser

    # The script usage
    usage  = "usage: %prog [options] experiment.json"

    # Create the parser
    self.parser = OptionParser(
      usage=usage,
      phil=self.phil_scope(),
      epilog=help_message)

    # Add an option to show configuration parameters
    self.parser.add_option(
      '-c',
      action='count',
      default=0,
      dest='show_config',
      help='Show the configuration parameters.')

  def phil_scope(self):
    ''' Get the phil scope. '''
    from libtbx.phil import parse
    new_phil_scope = parse('''
      integration {

        profile_model = 'profile_model.phil'
          .type = str
          .help = "The profile parameters output filename"

        integrated = 'integrated.pickle'
          .type = str
          .help = "The integrated output filename"

        reference = None
          .type = str
          .help = "The indexed reference spots input filename"

        predicted = None
          .type = str
          .help = "The predicted reflections input filename"
      }

      include scope dials.algorithms.integration.interface.phil_scope
      include scope dials.algorithms.profile_model.profile_model.phil_scope

    ''', process_includes=True)
    return new_phil_scope

  def run(self):
    ''' Perform the integration. '''
    from dials.util.command_line import heading
    from time import time

    # Check the number of arguments is correct
    start_time = time()

    # Parse the command line
    params, options, args = self.parser.parse_args()

    # Show config
    if options.show_config > 0:
      self.parser.print_phil(attributes_level=options.show_config-1)
      return

    # Check the number of command line arguments
    if len(args) != 1:
      self.parser.print_help()
      return

    print "=" * 80
    print ""
    print heading("Initialising")
    print ""

    # Print the diff phil
    diff_phil_str = self.parser.diff_phil().as_str()
    print 'Integrating with the following user specified parameters:\n'
    if (diff_phil_str is not ''):
      print diff_phil_str
    else:
      print 'All parameters set to defaults'

    # Load the experiment list
    experiments = self.load_experiments(args[0])

    # Load the data
    reference = self.load_reference(params.integration.reference)
    print ""
    predicted = self.load_predicted(params.integration.predicted)
    print ""

    # Initialise the integrator
    if None in experiments.goniometers():
      from dials.algorithms.integration import IntegratorStills
      integrator = IntegratorStills(params, experiments, reference, None, None)
    else:
      from dials.algorithms.profile_model.profile_model import ProfileModelList
      from dials.algorithms.integration.interface import IntegratorFactory
      from dials.array_family import flex

      # Compute the profile model
      # Predict the reflections
      # Match the predictions with the reference
      # Create the integrator
      if len(params.profile) > 1:
        assert(len(params.profile) == len(experiments) + 1)
        profile_model = ProfileModelList.load(params)
      else:
        assert(reference is not None)
        profile_model = ProfileModelList.compute(experiments, reference)
      print ""
      print "=" * 80
      print ""
      print heading("Predicting reflections")
      print ""
      if predicted is None:
        predicted = flex.reflection_table.from_predictions_multi(experiments)
      if reference:
        predicted.match_with_reference(reference)
      print ""
      integrator = IntegratorFactory.create(params, experiments, profile_model, predicted)

    # Integrate the reflections
    reflections = integrator.integrate()

    # Save the reflections
    self.save_reflections(reflections, params.integration.integrated)
    self.save_profile_model(profile_model, params.integration.profile_model)

    # Print the total time taken
    print "\nTotal time taken: ", time() - start_time

  def load_predicted(self, filename):
    ''' Load the predicted reflections. '''
    from dials.array_family import flex
    if filename is not None:
      return flex.reflection_table.from_pickle(filename)
    return None

  def load_experiments(self, filename):
    ''' Load the experiment list. '''
    from dxtbx.model.experiment.experiment_list import ExperimentListFactory
    from dials.util.command_line import Command
    Command.start('Loading experiments from %s' % filename)
    experiments = ExperimentListFactory.from_json_file(filename)
    Command.end('Loaded experiments from %s' % filename)
    if len(experiments) == 0:
      raise RuntimeError('experiment list is empty')
    elif len(experiments.imagesets()) > 1 or len(experiments.detectors()) > 1:
      raise RuntimeError('experiment list contains > 1 imageset or detector')
    return experiments

  def load_reference(self, filename):
    ''' Load the reference spots. '''
    from dials.util.command_line import Command
    from dials.array_family import flex

    if not filename:
      return None

    Command.start('Loading reference spots from %s' % filename)
    reference = flex.reflection_table.from_pickle(filename)
    assert("miller_index" in reference)
    Command.end('Loaded reference spots from %s' % filename)
    Command.start('Removing reference spots with invalid coordinates')
    mask = flex.bool([x == (0, 0, 0) for x in reference['xyzcal.mm']])
    reference.del_selected(mask)
    mask = flex.bool([h == (0, 0, 0) for h in reference['miller_index']])
    reference.del_selected(mask)
    Command.end('Removed reference spots with invalid coordinates, %d remaining' %
                len(reference))
    return reference

  def save_reflections(self, reflections, filename):
    ''' Save the reflections to file. '''
    from dials.util.command_line import Command
    Command.start('Saving %d reflections to %s' % (len(reflections), filename))
    reflections.as_pickle(filename)
    Command.end('Saved %d reflections to %s' % (len(reflections), filename))

  def save_profile_model(self, profile_model, filename):
    ''' Save the profile model parameters. '''
    from dials.util.command_line import Command
    Command.start('Saving the profile model parameters to %s' % filename)
    with open(filename, "w") as outfile:
      outfile.write(profile_model.dump().as_str())
    Command.end('Saved the profile model parameters to %s' % filename)


if __name__ == '__main__':
  script = Script()
  script.run()

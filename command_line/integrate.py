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

Examples::

  dials.integrate experiments.json reference=indexed.pickle

  dials.integrate experiments.json reference=indexed.pickle integrated=integrated.pickle

  dials.integrate experiments.json profile.phil

  dials.integrate experiments.json predicted=predicted.pickle reference=indexed.pickle

'''

# Create the phil scope
from libtbx.phil import parse
phil_scope = parse('''

  output {
    profile_model = 'profile_model.phil'
      .type = str
      .help = "The profile parameters output filename"

    reflections = 'integrated.pickle'
      .type = str
      .help = "The integrated output filename"
  }

  sampling
    .expert_level = 1
  {

    reflections_per_degree = 50
      .help = "The number of predicted reflections per degree of the sweep "
              "to integrate."
      .type = float(value_min=0.)

    minimum_sample_size = 1000
      .help = "cutoff that determines whether subsetting of the input "
              "prediction list is done"
      .type = int

    maximum_sample_size = None
      .help = "The maximum number of predictions to integrate."
              "Overrides reflections_per_degree if that produces a"
              "larger sample size."
      .type = int(value_min=1)

    integrate_all_reflections = True
      .help = "Override reflections_per_degree and integrate all predicted"
              "reflections."
      .type = bool

  }

  include scope dials.algorithms.integration.integrator.phil_scope
  include scope dials.algorithms.profile_model.factory.phil_scope
  include scope dials.algorithms.spot_prediction.reflection_predictor.phil_scope

''', process_includes=True)


class Script(object):
  ''' The integration program. '''

  def __init__(self):
    '''Initialise the script.'''
    from dials.util.options import OptionParser
    import libtbx.load_env

    # The script usage
    usage  = "usage: %s [options] experiment.json" % libtbx.env.dispatcher_name

    # Create the parser
    self.parser = OptionParser(
      usage=usage,
      phil=phil_scope,
      epilog=help_message,
      read_experiments=True,
      read_reflections=True)

  def run(self):
    ''' Perform the integration. '''
    from dials.util.command_line import heading
    from dials.util.options import flatten_reflections, flatten_experiments
    from time import time

    # Check the number of arguments is correct
    start_time = time()

    # Parse the command line
    params, options = self.parser.parse_args(show_diff_phil=True)
    reference = flatten_reflections(params.input.reflections)
    experiments = flatten_experiments(params.input.experiments)
    if len(reference) == 0:
      reference = None
    else:
      assert(len(reference) == 1)
      reference = reference[0]
    if len(experiments) == 0:
      self.parser.print_help()
      return
    elif len(experiments.imagesets()) > 1 or len(experiments.detectors()) > 1:
      raise RuntimeError('experiment list contains > 1 imageset or detector')

    print "=" * 80
    print ""
    print heading("Initialising")
    print ""

    # Load the data
    reference = self.process_reference(reference)
    print ""

    # Initialise the integrator
    # if None in experiments.goniometers():
    #   from dials.algorithms.integration import IntegratorStills
    #   integrator = IntegratorStills(params, experiments, reference, None, None)
    # else:
    if (True):
      from dials.algorithms.profile_model.factory import ProfileModelFactory
      from dials.algorithms.integration.integrator import IntegratorFactory
      from dials.array_family import flex

      # Compute the profile model
      # Predict the reflections
      # Match the predictions with the reference
      # Create the integrator
      profile_model = ProfileModelFactory.create(params, experiments, reference)
      print ""
      print "=" * 80
      print ""
      print heading("Predicting reflections")
      print ""
      predicted = profile_model.predict_reflections(
        experiments,
        dmin=params.prediction.dmin,
        dmax=params.prediction.dmax,
        margin=params.prediction.margin,
        force_static=params.prediction.force_static)

      if not params.sampling.integrate_all_reflections:
        nref_per_degree = params.sampling.reflections_per_degree
        min_sample_size = params.sampling.minimum_sample_size
        max_sample_size = params.sampling.maximum_sample_size

        # this code is very similar to David's code in algorithms/refinement/reflection_manager.py!

        # constants
        from math import pi
        RAD2DEG = 180. / pi
        DEG2RAD = pi / 180.

        working_isel = flex.size_t()
        for iexp, exp in enumerate(experiments):

          sel = predicted['id'] == iexp
          isel = sel.iselection()
          #refs = self._reflections.select(sel)
          nrefs = sample_size = len(isel)

          # set sample size according to nref_per_degree (per experiment)
          if exp.scan and nref_per_degree:
            sweep_range_rad = exp.scan.get_oscillation_range(deg=False)
            width = abs(sweep_range_rad[1] -
                        sweep_range_rad[0]) * RAD2DEG
            sample_size = int(nref_per_degree * width)
          else: sweep_range_rad = None

          # adjust sample size if below the chosen limit
          sample_size = max(sample_size, min_sample_size)

          # set maximum sample size if requested
          if max_sample_size:
            sample_size = min(sample_size, max_sample_size)

          # determine subset and collect indices
          if sample_size < nrefs:
            isel = isel.select(flex.random_selection(nrefs, sample_size))
          working_isel.extend(isel)

        # create subset
        predicted = predicted.select(working_isel)

      if reference:
        predicted.match_with_reference(reference)
      print ""
      integrator = IntegratorFactory.create(params, experiments, profile_model, predicted)

    # Integrate the reflections
    reflections = integrator.integrate()

    # Save the reflections
    self.save_reflections(reflections, params.output.reflections)
    self.save_profile_model(profile_model, params.output.profile_model)

    # Print the total time taken
    print "\nTotal time taken: ", time() - start_time

  def process_reference(self, reference):
    ''' Load the reference spots. '''
    from dials.util.command_line import Command
    from dials.array_family import flex
    if reference is None:
      return None
    assert("miller_index" in reference)
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
  from dials.util import halraiser
  try:
    script = Script()
    script.run()
  except Exception as e:
    halraiser(e)

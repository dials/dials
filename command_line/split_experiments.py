#!/usr/bin/env dials.python
from __future__ import absolute_import, division
from libtbx.phil import parse

help_message = """

Utility script to split experiments and reflections from single files into
multiple files with one experiment per output experiment file and one
reflection file per output experiment file.

Example::

  dials.split_experiments combined_experiments.json combined_reflections.pickle

"""

class Script(object):

  def __init__(self):
    '''Initialise the script.'''
    from dials.util.options import OptionParser
    from libtbx.phil import parse
    import libtbx.load_env

    # The phil scope
    phil_scope = parse('''
      by_detector = False
        .type = bool
        .help = "If True, instead of producing separate files for each"
                "experiment, experiments are grouped by unique detector"
                "model in the input set of experiments. For example, if"
                "there are five detector models in the input data, five"
                "sets of files will be produced, each containing"
                "experiments that reference a single detector model."

      output {
        experiments_prefix = experiments
          .type = str
          .help = "Filename prefix for the split experimental models"

        reflections_prefix = reflections
          .type = str
          .help = "Filename prefix for the split reflections"
      }
    ''', process_includes=True)

    # The script usage
    usage  = "usage: %s [options] [param.phil] " \
             "experiments1.json experiments2.json reflections1.pickle " \
             "reflections2.pickle..." \
             % libtbx.env.dispatcher_name

    # Create the parser
    self.parser = OptionParser(
      usage=usage,
      phil=phil_scope,
      read_reflections=True,
      read_experiments=True,
      check_format=False,
      epilog=help_message)

  def run(self):
    '''Execute the script.'''

    from dials.util.options import flatten_reflections, flatten_experiments
    from libtbx.utils import Sorry
    from dials.array_family import flex

    # Parse the command line
    params, options = self.parser.parse_args(show_diff_phil=True)

    # Try to load the models and data
    if not params.input.experiments:
      print "No Experiments found in the input"
      self.parser.print_help()
      return
    if params.input.reflections:
      if len(params.input.reflections) != len(params.input.experiments):
        raise Sorry("The number of input reflections files does not match the "
          "number of input experiments")

    experiments = flatten_experiments(params.input.experiments)
    if params.input.reflections:
      reflections = flatten_reflections(params.input.reflections)[0]
    else:
      reflections = None

    import math
    experiments_template = "%s_%%0%sd.json" %(
      params.output.experiments_prefix,
      int(math.floor(math.log10(len(experiments))) + 1))
    reflections_template = "%s_%%0%sd.pickle" %(
      params.output.reflections_prefix,
      int(math.floor(math.log10(len(experiments))) + 1))

    from dxtbx.model.experiment_list import ExperimentList
    from dxtbx.serialize import dump
    if params.by_detector:
      if reflections is None:
        split_data = {detector:{'experiments': ExperimentList()}
                      for detector in experiments.detectors()}
      else:
        split_data = {detector:{'experiments': ExperimentList(),
                                'reflections': flex.reflection_table()}
                      for detector in experiments.detectors()}

      for i, experiment in enumerate(experiments):
        split_expt_id = experiments.detectors().index(experiment.detector)
        experiment_filename = experiments_template % split_expt_id
        print 'Adding experiment %d to %s' %(i, experiment_filename)
        split_data[experiment.detector]['experiments'].append(experiment)
        if reflections is not None:
          reflections_filename = reflections_template % split_expt_id
          print 'Adding reflections for experiment %d to %s' %(i, reflections_filename)
          ref_sel = reflections.select(reflections['id'] == i)
          ref_sel['id'] = flex.int(len(ref_sel), len(split_data[experiment.detector]['experiments'])-1)
          split_data[experiment.detector]['reflections'].extend(ref_sel)

      for i, detector in enumerate(experiments.detectors()):
        experiment_filename = experiments_template %i
        print 'Saving experiment %d to %s' %(i, experiment_filename)
        dump.experiment_list(split_data[detector]['experiments'], experiment_filename)

        if reflections is not None:
          reflections_filename = reflections_template %i
          print 'Saving reflections for experiment %d to %s' %(i, reflections_filename)
          split_data[detector]['reflections'].as_pickle(reflections_filename)
    else:
      for i, experiment in enumerate(experiments):
        from dxtbx.model.experiment_list import ExperimentList
        from dxtbx.serialize import dump
        experiment_filename = experiments_template %i
        print 'Saving experiment %d to %s' %(i, experiment_filename)
        dump.experiment_list(ExperimentList([experiment]), experiment_filename)

        if reflections is not None:
          reflections_filename = reflections_template %i
          print 'Saving reflections for experiment %d to %s' %(i, reflections_filename)
          ref_sel = reflections.select(reflections['id'] == i)
          ref_sel['id'] = flex.int(len(ref_sel), 0)
          ref_sel.as_pickle(reflections_filename)

    return

if __name__ == "__main__":
  from dials.util import halraiser
  try:
    script = Script()
    script.run()
  except Exception as e:
    halraiser(e)

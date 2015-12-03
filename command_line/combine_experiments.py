#!/usr/bin/env dials.python
from __future__ import division
import sys, os
from libtbx.phil import command_line, parse

help_message = '''

Utility script to combine multiple reflections and experiments files into
one multi-experiment reflections and one experiments file. Experiments are
matched to reflections in the order they are provided as input.

Reference models can be chosen from any of the input experiments files. These
will replace all other models of that type in the output experiments file.
This is useful, for example, for combining mutiple experiments that should
differ only in their crystal models. No checks are made to ensure that a
reference model is a good replacement model.

Although only one reference model of each type is allowed, more complex
combinations of experiments can be created by repeat runs.

Examples::

  dials.combine_experiments experiments_0.json experiments_1.json \\
    reflections_0.pickle reflections_1.pickle \\
    reference_from_experiment.beam=0 \\
    reference_from_experiment.detector=0

'''

class CombineWithReference(object):

  def __init__(self, beam=None, goniometer=None, scan=None,
                     crystal=None, detector=None):

    self.ref_beam = beam
    self.ref_goniometer = goniometer
    self.ref_scan = scan
    self.ref_crystal = crystal
    self.ref_detector = detector

    return

  def __call__(self, experiment):

    beam = experiment.beam if self.ref_beam is None else self.ref_beam
    goniometer = experiment.goniometer if self.ref_goniometer is None \
      else self.ref_goniometer
    scan = experiment.scan if self.ref_scan is None else self.ref_scan
    crystal = experiment.crystal if self.ref_crystal is None \
      else self.ref_crystal
    detector = experiment.detector if self.ref_detector is None \
      else self.ref_detector

    from dxtbx.model.experiment.experiment_list import Experiment
    return Experiment(beam=beam,
                      detector=detector,
                      scan=scan,
                      goniometer=goniometer,
                      crystal=crystal,
                      imageset=experiment.imageset)

class Script(object):

  def __init__(self):
    '''Initialise the script.'''
    from dials.util.options import OptionParser
    from libtbx.phil import parse
    import libtbx.load_env

    # The phil scope
    phil_scope = parse('''

      reference_from_experiment{
        beam = None
          .help = "Take beam model from this experiment to overwrite all other"
                  "beam models in the combined experiments"
          .type = int(value_min=0)

        scan = None
          .help = "Take scan model from this experiment to overwrite all other"
                  "scan models in the combined experiments"
          .type = int(value_min=0)

        crystal = None
          .help = "Take crystal model from this experiment to overwrite all"
                  "other crystal models in the combined experiments"
          .type = int(value_min=0)

        goniometer = None
          .help = "Take goniometer model from this experiment to overwrite all"
                  "other goniometer models in the combined experiments"
          .type = int(value_min=0)

        detector = None
          .help = "Take detector model from this experiment to overwrite all"
                  "other detector models in the combined experiments"
          .type = int(value_min=0)
      }

      output {
        experiments_filename = combined_experiments.json
          .type = str
          .help = "The filename for combined experimental models"

        reflections_filename = combined_reflections.pickle
          .type = str
          .help = "The filename for combined reflections"
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

    # Parse the command line
    params, options = self.parser.parse_args(show_diff_phil=True)

    # Try to load the models and data
    if len(params.input.experiments) == 0:
      print "No Experiments found in the input"
      self.parser.print_help()
      return
    if len(params.input.reflections) == 0:
      print "No reflection data found in the input"
      self.parser.print_help()
      return
    try:
      assert len(params.input.reflections) == len(params.input.experiments)
    except AssertionError:
      raise Sorry("The number of input reflections files does not match the "
        "number of input experiments")

    flat_exps = flatten_experiments(params.input.experiments)

    ref_beam = params.reference_from_experiment.beam
    ref_goniometer = params.reference_from_experiment.goniometer
    ref_scan = params.reference_from_experiment.scan
    ref_crystal = params.reference_from_experiment.crystal
    ref_detector = params.reference_from_experiment.detector

    if ref_beam is not None:
      try:
        ref_beam = flat_exps[ref_beam].beam
      except IndexError:
        raise Sorry("{0} is not a valid experiment ID".format(ref_beam))

    if ref_goniometer is not None:
      try:
        ref_goniometer = flat_exps[ref_goniometer].goniometer
      except IndexError:
        raise Sorry("{0} is not a valid experiment ID".format(ref_goniometer))

    if ref_scan is not None:
      try:
        ref_scan = flat_exps[ref_scan].scan
      except IndexError:
        raise Sorry("{0} is not a valid experiment ID".format(ref_scan))

    if ref_crystal is not None:
      try:
        ref_crystal = flat_exps[ref_crystal].crystal
      except IndexError:
        raise Sorry("{0} is not a valid experiment ID".format(ref_crystal))

    if ref_detector is not None:
      try:
        ref_detector = flat_exps[ref_detector].detector
      except IndexError:
        raise Sorry("{0} is not a valid experiment ID".format(ref_detector))

    combine = CombineWithReference(beam=ref_beam, goniometer=ref_goniometer,
                  scan=ref_scan, crystal=ref_crystal, detector=ref_detector)

    # set up global experiments and reflections lists
    from dials.array_family import flex
    reflections = flex.reflection_table()
    global_id = 0
    from dxtbx.model.experiment.experiment_list import ExperimentList
    experiments=ExperimentList()

    # loop through the input, building up the global lists
    nrefs_per_exp = []
    for ref_wrapper, exp_wrapper in zip(params.input.reflections,
                                        params.input.experiments):
      refs = ref_wrapper.data
      exps = exp_wrapper.data
      for i, exp in enumerate(exps):
        sel = refs['id'] == i
        sub_ref = refs.select(sel)
        nrefs_per_exp.append(len(sub_ref))
        sub_ref['id'] = flex.int(len(sub_ref), global_id)
        reflections.extend(sub_ref)
        experiments.append(combine(exp))
        global_id += 1

    # print number of reflections per experiment
    from libtbx.table_utils import simple_table
    header = ["Experiment", "Nref"]
    rows = [(str(i), str(n)) for (i, n) in enumerate(nrefs_per_exp)]
    st = simple_table(rows, header)
    print st.format()

    # save output
    from dxtbx.model.experiment.experiment_list import ExperimentListDumper
    print 'Saving combined experiments to {0}'.format(
      params.output.experiments_filename)
    dump = ExperimentListDumper(experiments)
    dump.as_json(params.output.experiments_filename)
    print 'Saving combined reflections to {0}'.format(
      params.output.reflections_filename)
    reflections.as_pickle(params.output.reflections_filename)

    return

if __name__ == "__main__":
  from dials.util import halraiser
  try:
    script = Script()
    script.run()
  except Exception as e:
    halraiser(e)


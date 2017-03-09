#!/usr/bin/env dials.python
from __future__ import absolute_import, division
from libtbx.phil import parse

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

    average_detector = False
      .help = "Create an average detector model from all the input detector"
              "models and use it as the reference. Not compatible with"
              "reference_from_experiment.detector"
      .type = bool

    compare_models = True
      .help = "Whether to compare a model with the reference model before"
              "replacing it. If the comparison falls outside the tolerance,"
              "the combination will not be allowed. Disable comparison to force"
              "overwriting of models with the reference"
      .type = bool

    average_hierarchy_level = None
      .help = "For hierarchical detectors, optionally provide a single level"
              "to do averaging at."
      .type = int(value_min=0)

    include scope dials.util.options.tolerance_phil_scope

  }

  output {
    experiments_filename = combined_experiments.json
      .type = str
      .help = "The filename for combined experimental models"

    reflections_filename = combined_reflections.pickle
      .type = str
      .help = "The filename for combined reflections"

    n_subset = None
      .type = int
      .help = "If not None, keep a random subset of size n_subset when"
              "saving the combined experiments"

    max_batch_size = None
      .type = int
      .expert_level = 2
      .help = "If not None, split the resultant combined set of experiments"
              "into seperate files, each at most max_batch_size number of"
              "experiments. Example, if there were 5500 experiments and"
              "max_batch_size is 1000, 6 experiment lists will be created,"
              "of sizes 917, 917, 917, 917, 916, 916"

    delete_shoeboxes = False
      .type = bool
      .expert_level = 2
      .help = "If true, delete shoeboxes from reflection tables while comb-"
              "ining them to save on memory."
  }
''', process_includes=True)

class CombineWithReference(object):

  def __init__(self, beam=None, goniometer=None, scan=None,
                     crystal=None, detector=None, params=None):

    self.ref_beam = beam
    self.ref_goniometer = goniometer
    self.ref_scan = scan
    self.ref_crystal = crystal
    self.ref_detector = detector
    self.tolerance = None
    if params:
      if params.reference_from_experiment.compare_models:
        self.tolerance = params.reference_from_experiment.tolerance
      self.average_detector = params.reference_from_experiment.average_detector
    else:
      self.average_detector = False

    return

  def __call__(self, experiment):
    from dxtbx.datablock import BeamComparison
    from dxtbx.datablock import DetectorComparison
    from dxtbx.datablock import GoniometerComparison

    if self.tolerance:
      compare_beam = BeamComparison(
        wavelength_tolerance=self.tolerance.beam.wavelength,
        direction_tolerance=self.tolerance.beam.direction,
        polarization_normal_tolerance=self.tolerance.beam.polarization_normal,
        polarization_fraction_tolerance=self.tolerance.beam.polarization_fraction)
      compare_detector = DetectorComparison(
        fast_axis_tolerance=self.tolerance.detector.fast_axis,
        slow_axis_tolerance=self.tolerance.detector.slow_axis,
        origin_tolerance=self.tolerance.detector.origin)
      compare_goniometer = GoniometerComparison(
        rotation_axis_tolerance=self.tolerance.goniometer.rotation_axis,
        fixed_rotation_tolerance=self.tolerance.goniometer.fixed_rotation,
        setting_rotation_tolerance=self.tolerance.goniometer.setting_rotation)

    else:
      compare_beam = None
      compare_detector = None
      compare_goniometer = None

    if self.ref_beam:
      if compare_beam:
        assert(compare_beam(self.ref_beam, experiment.beam))
      beam = self.ref_beam
    else:
      beam = experiment.beam

    if self.ref_detector and self.average_detector:
      detector = self.ref_detector
    elif self.ref_detector and not self.average_detector:
      if compare_detector:
        assert(compare_detector(self.ref_detector, experiment.detector))
      detector = self.ref_detector
    else:
      detector = experiment.detector

    if self.ref_goniometer:
      if compare_goniometer:
        assert(compare_goniometer(self.ref_goniometer, experiment.goniometer))
      goniometer = self.ref_goniometer
    else:
      goniometer = experiment.goniometer

    if self.ref_scan:
      scan = self.ref_scan
    else:
      scan = experiment.scan

    if self.ref_crystal:
      crystal = self.ref_crystal
    else:
      crystal = experiment.crystal

    from dxtbx.model.experiment_list import Experiment
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
    import libtbx.load_env

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

    from dials.util.options import flatten_experiments
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
      assert not params.reference_from_experiment.average_detector
      try:
        ref_detector = flat_exps[ref_detector].detector
      except IndexError:
        raise Sorry("{0} is not a valid experiment ID".format(ref_detector))
    elif params.reference_from_experiment.average_detector:
      # Average all of the detectors together
      from scitbx.matrix import col
      def average_detectors(target, panelgroups, depth):
        # Recursive function to do the averaging

        if params.reference_from_experiment.average_hierarchy_level is None or \
            depth == params.reference_from_experiment.average_hierarchy_level:
          n = len(panelgroups)
          sum_fast = col((0.0,0.0,0.0))
          sum_slow = col((0.0,0.0,0.0))
          sum_ori  = col((0.0,0.0,0.0))

          # Average the d matrix vectors
          for pg in panelgroups:
            sum_fast += col(pg.get_local_fast_axis())
            sum_slow += col(pg.get_local_slow_axis())
            sum_ori  += col(pg.get_local_origin())
          sum_fast /= n
          sum_slow /= n
          sum_ori  /= n

          # Re-orthagonalize the slow and the fast vectors by rotating around the cross product
          c = sum_fast.cross(sum_slow)
          a = sum_fast.angle(sum_slow, deg=True)/2
          sum_fast = sum_fast.rotate(c, a-45, deg=True)
          sum_slow = sum_slow.rotate(c, -(a-45), deg=True)

          target.set_local_frame(sum_fast,sum_slow,sum_ori)

        if target.is_group():
          # Recurse
          for i, target_pg in enumerate(target):
            average_detectors(target_pg, [pg[i] for pg in panelgroups], depth+1)

      ref_detector = flat_exps[0].detector
      average_detectors(ref_detector.hierarchy(), [e.detector.hierarchy() for e in flat_exps], 0)

    combine = CombineWithReference(beam=ref_beam, goniometer=ref_goniometer,
                  scan=ref_scan, crystal=ref_crystal, detector=ref_detector,
                  params=params)

    # set up global experiments and reflections lists
    from dials.array_family import flex
    reflections = flex.reflection_table()
    global_id = 0
    from dxtbx.model.experiment_list import ExperimentList
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
        if params.output.delete_shoeboxes and 'shoebox' in sub_ref:
          del sub_ref['shoebox']
        reflections.extend(sub_ref)
        experiments.append(combine(exp))
        global_id += 1

    # print number of reflections per experiment
    from libtbx.table_utils import simple_table
    header = ["Experiment", "Nref"]
    rows = [(str(i), str(n)) for (i, n) in enumerate(nrefs_per_exp)]
    st = simple_table(rows, header)
    print st.format()

    # save a random subset if requested
    if params.output.n_subset is not None and len(experiments) > params.output.n_subset:
      import random
      subset_exp = ExperimentList()
      subset_refls = flex.reflection_table()
      n_picked = 0
      indices = range(len(experiments))
      while n_picked < params.output.n_subset:
        idx = indices.pop(random.randint(0, len(indices)-1))
        subset_exp.append(experiments[idx])
        refls = reflections.select(reflections['id'] == idx)
        refls['id'] = flex.int(len(refls), n_picked)
        subset_refls.extend(refls)
        n_picked += 1
      experiments = subset_exp
      reflections = subset_refls

    def save_output(experiments, reflections, exp_name, refl_name):
      # save output
      from dxtbx.model.experiment_list import ExperimentListDumper
      print 'Saving combined experiments to {0}'.format(exp_name)
      dump = ExperimentListDumper(experiments)
      dump.as_json(exp_name)
      print 'Saving combined reflections to {0}'.format(refl_name)
      reflections.as_pickle(refl_name)

    if params.output.max_batch_size is None:
      save_output(experiments, reflections, params.output.experiments_filename, params.output.reflections_filename)
    else:
      from dxtbx.command_line.image_average import splitit
      import os
      for i, indices in enumerate(splitit(range(len(experiments)), (len(experiments)//params.output.max_batch_size)+1)):
        batch_expts = ExperimentList()
        batch_refls = flex.reflection_table()
        for sub_id, sub_idx in enumerate(indices):
          batch_expts.append(experiments[sub_idx])
          sub_refls = reflections.select(reflections['id'] == sub_idx)
          sub_refls['id'] = flex.int(len(sub_refls), sub_id)
          batch_refls.extend(sub_refls)
        exp_filename = os.path.splitext(params.output.experiments_filename)[0] + "_%03d.json"%i
        ref_filename = os.path.splitext(params.output.reflections_filename)[0] + "_%03d.pickle"%i
        save_output(batch_expts, batch_refls, exp_filename, ref_filename)
    return

if __name__ == "__main__":
  from dials.util import halraiser
  try:
    script = Script()
    script.run()
  except Exception as e:
    halraiser(e)

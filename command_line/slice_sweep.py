#!/usr/bin/env python
#
# dials.slice_sweep.py
#
#  Copyright (C) 2014 Diamond Light Source and STFC Rutherford Appleton
#  Laboratory, UK.
#
#  Author: David Waterman
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division
from os.path import splitext, basename
from scitbx.array_family import flex
from libtbx.utils import Sorry

help_message = '''

Slice a sweep to produce a smaller sweep within the bounds of the original. If
experiments are provided, modify the scan objects within the experiments. If
reflections are provided, remove reflections outside the provided scan ranges.
Each scan_range parameter refers to a single experiment ID, counting up from
zero. Any reflections with experiment ID not matched by a scan_range parameter
are are removed.

Examples::

  dials.slice_sweep experiments.json reflections.pickle "scan_range=1 20"

  # two experiments and reflections with IDs '0' and '1'
  dials.slice_sweep experiments.json reflections.pickle \
    "scan_range=1 20" "scan_range=5 30"

'''

# The phil scope
from libtbx.phil import parse
phil_scope = parse('''

  output {
    experiments_filename = None
      .type = str
      .help = "The filename for output experimental models with sliced scans.
               By default generated automatically from the input name"

    reflections_filename = None
      .type = str
      .help = "The filename for output reflections sliced to remove those"
              "outside the reduced scan range. By default generated"
              "automatically from the input name"
  }

  scan_range = None
    .type = str
    .help = "Scan range in images to slice a sweep. Number of arguments"
            "must be a factor of two. Specifying \"0 0\" will use all images"
            "by default. The given range follows C conventions"
            "(e.g. j0 <= j < j1)."
    .type = ints(size=2)
    .multiple = True

''')
#''', process_includes=True)

def slice_experiments(experiments, scan_ranges):
  '''

  :param experiments
  :type experiments: dxtbx.model.experiment.experiment_list.ExperimentList
  :param scan_range:
  :type scan_range: list of 2-tuples defining scan range for each experiment'''

  # copy the experiments
  import copy
  experiments = copy.deepcopy(experiments)

  if len(experiments) != len(scan_ranges):
    raise Sorry("Input experiment list and scan_ranges are not of the same length")

  for exp, sr in zip(experiments, scan_ranges):
    if sr is None: continue
    im_range = exp.scan.get_image_range()
    if sr[0] < im_range[0] or sr[1] > im_range[1]:
      raise IndexError("requested slice outside current scan range")

    exp.scan.set_image_range(sr)

  return experiments

def slice_reflections(reflections, scan_ranges):
  '''

  :param reflections: reflection table of input reflections which must contain
                      the
  :type reflections: dials.array_family.flex.reflection_table
  :param scan_range: list of 2-tuples defining scan range for each experiment
                     id contained within the reflections
  :type scan_range: list of 2-tuples defining scan range for each experiment'''

  # copy the reflections
  import copy
  reflections = copy.deepcopy(reflections)

  to_keep = flex.size_t()
  for iexp, sr in enumerate(scan_ranges):

    if sr is None: continue
    isel = (reflections['id'] == iexp).iselection()
    frames = (reflections['xyzobs.px.value'].parts()[2]).select(isel)
    # reflns on image n have frames in range [n-1, n)
    in_low_lim = frames >= sr[0] - 1
    in_high_lim = frames < sr[1]
    in_lim = in_low_lim & in_high_lim

    # which indices to keep?
    sub_isel = isel.select(in_lim)
    to_keep.extend(sub_isel)

  # implictly also removes any reflections with ID outside the range of the
  # length of scan_ranges
  return reflections.select(to_keep)

def calculate_frame_numbers(reflections, experiments):
  """calculate observed frame numbers for all reflections, if not already
  set"""

  # FIXME this is adapted from ReflectionManager code. It does not really
  # need to exist twice. Could ReflectionManager's version become a
  # staticmethod instead, and then just reuse that?

  # Only do this if we have to
  if reflections.has_key('xyzobs.px.value'): return

  # Ok, frames are not set, so set them, with dummy observed pixel values
  frames = flex.double(len(reflections), 0.)
  for iexp, exp in enumerate(experiments):
    scan = exp.scan
    if not scan: continue
    sel = reflections['id'] == iexp
    xyzobs = reflections["xyzobs.mm.value"].select(sel)
    angles = xyzobs.parts()[2]
    to_update = scan.get_array_index_from_angle(angles, deg=False)
    frames.set_selected(sel, to_update)
  reflections['xyzobs.px.value'] = flex.vec3_double(
          flex.double(len(reflections), 0.),
          flex.double(len(reflections), 0.),
          frames)

  return reflections

class Script(object):
  '''A class for running the script.'''

  def __init__(self):
    '''Initialise the script.'''
    from dials.util.options import OptionParser
    import libtbx.load_env

    # The script usage
    usage  = "usage: %s [options] [param.phil] " \
             "experiments.json reflections.pickle" \
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
    import cPickle as pickle

    # Parse the command line
    params, options = self.parser.parse_args(show_diff_phil=True)
    reflections = flatten_reflections(params.input.reflections)
    experiments = flatten_experiments(params.input.experiments)

    # Try to load the models and data
    slice_exps = len(experiments) > 0
    slice_refs = len(reflections) > 0

    if len(reflections) > 1:
      raise Sorry("Only one reflections list can be imported at present")
    reflections = reflections[0]

    # calculate frame numbers if needed
    if len(experiments) > 0:
      reflections = calculate_frame_numbers(reflections, experiments)

    # if we still don't have the right column give up
    if not reflections.has_key('xyzobs.px.value'):
      raise Sorry("These reflections do not have frame numbers set, and "
        "there are no experiments provided to calculate these.")

    # set trivial case where no scan range is provided at all
    if len(params.scan_range) == 0:
      params.scan_range = [None]

    # do slicing
    if slice_exps:
      sliced_experiments = slice_experiments(experiments, params.scan_range)
    if slice_refs:
      sliced_reflections = slice_reflections(reflections, params.scan_range)

    # Save sliced experiments
    if slice_exps:
      output_experiments_filename = params.output.experiments_filename
      if output_experiments_filename is None:
        # take first filename as template
        bname = basename(params.input.experiments[0].filename)
        bname = splitext(bname)[0]
        if not bname: bname = "experiments"
        if len(params.scan_range) == 1 and params.scan_range[0] is not None:
          ext = "_{0}_{1}.json".format(*params.scan_range[0])
        else:
          ext = "_subsets.json"
        output_experiments_filename = bname + ext
      print 'Saving sliced experiments to {0}'.format(
        output_experiments_filename)

      from dxtbx.model.experiment.experiment_list import ExperimentListDumper
      dump = ExperimentListDumper(sliced_experiments)
      dump.as_json(output_experiments_filename)

    # Save sliced reflections
    if slice_refs:
      output_reflections_filename = params.output.reflections_filename
      if output_reflections_filename is None:
        # take first filename as template
        bname = basename(params.input.reflections[0].filename)
        bname = splitext(bname)[0]
        if not bname: bname = "reflections"
        if len(params.scan_range) == 1 and params.scan_range[0] is not None:
          ext = "_{0}_{1}.pickle".format(*params.scan_range[0])
        else:
          ext = "_subsets.pickle"
        output_reflections_filename = bname + ext

      print 'Saving sliced reflections to {0}'.format(
        output_reflections_filename)
      sliced_reflections.as_pickle(output_reflections_filename)

    return

if __name__ == '__main__':
  from dials.util import halraiser
  try:
    script = Script()
    script.run()
  except Exception as e:
    halraiser(e)

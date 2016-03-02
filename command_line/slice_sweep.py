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
from dials.algorithms.refinement.refinement_helpers import \
  calculate_frame_numbers

help_message = '''

Slice a sweep to produce a smaller sweep within the bounds of the original. If
experiments or datablocks are provided, modify the scan objects within these. If
reflections are provided, remove reflections outside the provided scan ranges.
Each scan_range parameter refers to a single experiment ID, counting up from
zero. Any reflections with experiment ID not matched by a scan_range parameter
are removed.

Examples::

  dials.slice_sweep experiments.json reflections.pickle "scan_range=1 20"

  dials.slice_sweep datablock.json "scan_range=1 20"

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

    datablocks_filename = None
      .type = str
      .help = "The filename for the output datablock with sliced scans.
               By default generated automatically from the input name"

  }

  scan_range = None
    .type = str
    .help = "Scan range in images to slice a sweep. Number of arguments"
            "must be a factor of two. Specifying \"0 0\" will use all images"
            "by default. The given range follows C conventions"
            "(e.g. j0 <= j < j1)."
    .type = ints(size=2)
    .multiple = True

  block_size = None
    .type = float
    .help = "Overrides scan_range if present. This option splits each sweep"
            "into the nearest integer number of equal size blocks close to"
            "block_size degrees in width"

''')

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

    # slicing uses the array range, not the image range
    beg = sr[0] - 1
    end = sr[1]
    exp.scan.swap(exp.scan[beg:end])

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

def slice_datablocks(datablocks, scan_ranges):
  '''

  :param datablocks:
  :type datablocks: list of dxtbx.datablock.DataBlocks
  :param scan_range:
  :type scan_range: list of 2-tuples defining scan range for each unique scan'''

  # copy the experiments
  import copy
  datablocks = copy.deepcopy(datablocks)

  scans = []
  for db in datablocks:
    scans.extend(db.unique_scans())

  if len(scans) != len(scan_ranges):
    raise Sorry("Unique scans in the input datablock list and supplied scan_ranges are not of the same length")

  for scan, sr in zip(scans, scan_ranges):
    if sr is None: continue

    im_range = scan.get_image_range()
    if sr[0] < im_range[0] or sr[1] > im_range[1]:
      raise IndexError("requested slice outside current scan range")

    # slicing uses the array range, not the image range
    beg = sr[0] - 1
    end = sr[1]
    scan.swap(scan[beg:end])

  return datablocks

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
      read_datablocks=True,
      check_format=False,
      epilog=help_message)

  def run(self):
    '''Execute the script.'''

    from dials.util.options import flatten_reflections, flatten_experiments, \
      flatten_datablocks
    import cPickle as pickle

    # Parse the command line
    params, options = self.parser.parse_args(show_diff_phil=True)
    reflections = flatten_reflections(params.input.reflections)
    experiments = flatten_experiments(params.input.experiments)
    datablocks = flatten_datablocks(params.input.datablock)

    # Try to load the models and data
    slice_exps = len(experiments) > 0
    slice_refs = len(reflections) > 0
    slice_dbs = len(datablocks) > 0

    # Catch case of nothing to do
    if not any([slice_exps, slice_refs, slice_dbs]):
      print "No suitable input provided"
      self.parser.print_help()
      return

    if reflections:
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
    if slice_dbs:
      sliced_datablocks =  slice_datablocks(datablocks, params.scan_range)

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

    # Save sliced datablocks
    if slice_dbs:
      output_datablocks_filename = params.output.datablocks_filename
      if output_datablocks_filename is None:
        # take first filename as template
        bname = basename(params.input.datablock[0].filename)
        bname = splitext(bname)[0]
        if not bname: bname = "datablock"
        if len(params.scan_range) == 1 and params.scan_range[0] is not None:
          ext = "_{0}_{1}.json".format(*params.scan_range[0])
        else:
          ext = "_subsets.json"
        output_datablocks_filename = bname + ext
      print 'Saving sliced datablocks to {0}'.format(
        output_datablocks_filename)

      from dxtbx.datablock import DataBlockDumper
      dump = DataBlockDumper(sliced_datablocks)
      dump.as_file(output_datablocks_filename)

    return

if __name__ == '__main__':
  from dials.util import halraiser
  try:
    script = Script()
    script.run()
  except Exception as e:
    halraiser(e)

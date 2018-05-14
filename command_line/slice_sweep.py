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

from __future__ import absolute_import, division, print_function

from os.path import basename, splitext

from dials.algorithms.refinement.refinement_helpers import \
    calculate_frame_numbers
from dxtbx.datablock import DataBlock
from dxtbx.model.experiment_list import ExperimentList
from libtbx.utils import Sorry
from scitbx.array_family import flex

help_message = '''

Slice a sweep to produce a smaller sweep within the bounds of the original. If
experiments or datablocks are provided, modify the scan objects within these. If
reflections are provided, remove reflections outside the provided image ranges.
Each image_range parameter refers to a single experiment ID, counting up from
zero. Any reflections with experiment ID not matched by a image_range parameter
are removed.

Examples::

  dials.slice_sweep experiments.json reflections.pickle "image_range=1 20"

  dials.slice_sweep datablock.json "image_range=1 20"

  # two experiments and reflections with IDs '0' and '1'
  dials.slice_sweep experiments.json reflections.pickle \
    "image_range=1 20" "image_range=5 30"

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
              "outside the reduced image range. By default generated"
              "automatically from the input name"

    datablocks_filename = None
      .type = str
      .help = "The filename for the output datablock with sliced scans.
               By default generated automatically from the input name"

  }

  image_range = None
    .help = "Range in images to slice a sweep. The number of arguments"
            "must be a factor of two. Each pair of arguments gives a range"
            "that follows C conventions (e.g. j0 <= j < j1) when slicing the"
            "reflections by observed centroid."
    .type = ints(size=2)
    .multiple = True

  block_size = None
    .type = float
    .help = "Overrides image_range if present. This option splits each sweep"
            "into the nearest integer number of equal size blocks close to"
            "block_size degrees in width"

''')


def calculate_block_ranges(scan, block_size):
  '''

  :param scans
  :type a scan object
  :param block_size:
  :type block_size: target block size in degrees'''

  image_ranges = []
  nimages = scan.get_num_images()
  osc_range = scan.get_oscillation_range(deg=True)
  osc_width = abs(osc_range[1] - osc_range[0])
  nblocks = max(int(round(osc_width / block_size)), 1)
  nblocks = min(nblocks, nimages)
  # equal sized blocks except the last one that may contain extra images
  # to make up the remainder
  nimages_per_block = [nimages // nblocks] * (nblocks - 1) + \
    [nimages // nblocks + nimages % nblocks]
  start = scan.get_image_range()[0]
  for nim in nimages_per_block:
    image_ranges.append((start, start + nim - 1))
    start += nim

  return image_ranges

def slice_experiments(experiments, image_ranges):
  '''

  :param experiments
  :type experiments: dxtbx.model.experiment_list.ExperimentList
  :param image_range:
  :type image_range: list of 2-tuples defining scan range for each experiment'''

  # copy the experiments
  import copy
  experiments = copy.deepcopy(experiments)

  if len(experiments) != len(image_ranges):
    raise Sorry("Input experiment list and image_ranges are not of the same length")

  for exp, sr in zip(experiments, image_ranges):
    if sr is None: continue
    im_range = exp.scan.get_image_range()
    if sr[0] < im_range[0] or sr[1] > im_range[1]:
      raise IndexError("requested slice outside current scan range")

    # slicing uses the array range, not the image range
    arr_start = exp.scan.get_array_range()[0]
    beg = sr[0] - 1 - arr_start
    end = sr[1] - arr_start
    exp.scan.swap(exp.scan[beg:end])

  return experiments

def slice_reflections(reflections, image_ranges):
  '''

  :param reflections: reflection table of input reflections
  :type reflections: dials.array_family.flex.reflection_table
  :param image_range: list of 2-tuples defining scan range for each experiment
                     id contained within the reflections
  :type image_range: list of 2-tuples defining scan range for each experiment'''

  # copy the reflections
  import copy
  reflections = copy.deepcopy(reflections)

  to_keep = flex.size_t()
  for iexp, sr in enumerate(image_ranges):

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
  # length of image_ranges
  return reflections.select(to_keep)

def slice_datablocks(datablocks, image_ranges):
  '''

  :param datablocks:
  :type datablocks: list of dxtbx.datablock.DataBlocks
  :param image_range:
  :type image_range: list of 2-tuples defining scan range for each unique scan'''

  # copy the experiments
  import copy
  datablocks = copy.deepcopy(datablocks)

  scans = []
  for db in datablocks:
    scans.extend(db.unique_scans())

  if len(scans) != len(image_ranges):
    raise Sorry("Unique scans in the input datablock list and supplied image_ranges are not of the same length")

  for scan, sr in zip(scans, image_ranges):
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
      print("No suitable input provided")
      self.parser.print_help()
      return

    if reflections:
      if len(reflections) > 1:
        raise Sorry("Only one reflections list can be imported at present")
      reflections = reflections[0]

      # calculate frame numbers if needed
      if experiments:
        reflections = calculate_frame_numbers(reflections, experiments)

      # if we still don't have the right column give up
      if 'xyzobs.px.value' not in reflections:
        raise Sorry("These reflections do not have frame numbers set, and "
          "there are no experiments provided to calculate these.")

    # set trivial case where no scan range is provided at all
    if not params.image_range:
      params.image_range = [None]

    # check if slicing into blocks
    if params.block_size is not None:
      # in this case for simplicity, ensure that there is either an
      # an experiment list or datablocks, but not both. Ensure there is only
      # a single scan contained within.
      if [slice_exps, slice_dbs].count(True) != 1:
        raise Sorry("For slicing into blocks please provide either datablocks"
          " or experiments, but not both.")
      if slice_exps:
        if len(experiments) > 1:
          raise Sorry("For slicing into blocks please provide a single "
                      "scan only")
        scan = experiments[0].scan
      if slice_dbs:
        scans = datablocks[0].unique_scans()
        if len(scans) > 1 or len(datablocks) > 1:
          raise Sorry("For slicing into blocks please provide a single "
                      "scan only")
        scan = scans[0]

      # Having extracted the scan, calculate the blocks
      params.image_range = calculate_block_ranges(scan, params.block_size)

      # Do the slicing then recombine
      if slice_exps:
        sliced = [slice_experiments(experiments, [sr])[0] \
          for sr in params.image_range]
        sliced_experiments = ExperimentList()
        for exp in sliced:
          sliced_experiments.append(exp)

      if slice_dbs:
        sliced = [slice_datablocks(datablocks, [sr])[0] \
          for sr in params.image_range]
        imagesets = [db.extract_imagesets()[0] for db in sliced]
        sliced_datablocks = DataBlock(imagesets)

      # slice reflections if present
      if slice_refs:
        sliced = [slice_reflections(reflections, [sr]) \
          for sr in params.image_range]
        sliced_reflections = sliced[0]
        for i, rt in enumerate(sliced[1:]):
          rt['id'] += (i + 1) # set id
          sliced_reflections.extend(rt)

    else:
      # slice each dataset into the requested subset
      if slice_exps:
        sliced_experiments = slice_experiments(experiments, params.image_range)
      if slice_refs:
        sliced_reflections = slice_reflections(reflections, params.image_range)
      if slice_dbs:
        sliced_datablocks = slice_datablocks(datablocks, params.image_range)

    # Save sliced experiments
    if slice_exps:
      output_experiments_filename = params.output.experiments_filename
      if output_experiments_filename is None:
        # take first filename as template
        bname = basename(params.input.experiments[0].filename)
        bname = splitext(bname)[0]
        if not bname: bname = "experiments"
        if len(params.image_range) == 1 and params.image_range[0] is not None:
          ext = "_{0}_{1}.json".format(*params.image_range[0])
        else:
          ext = "_sliced.json"
        output_experiments_filename = bname + ext
      print('Saving sliced experiments to {0}'.format(
        output_experiments_filename))

      from dxtbx.model.experiment_list import ExperimentListDumper
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
        if len(params.image_range) == 1 and params.image_range[0] is not None:
          ext = "_{0}_{1}.pickle".format(*params.image_range[0])
        else:
          ext = "_sliced.pickle"
        output_reflections_filename = bname + ext

      print('Saving sliced reflections to {0}'.format(
        output_reflections_filename))
      sliced_reflections.as_pickle(output_reflections_filename)

    # Save sliced datablocks
    if slice_dbs:
      output_datablocks_filename = params.output.datablocks_filename
      if output_datablocks_filename is None:
        # take first filename as template
        bname = basename(params.input.datablock[0].filename)
        bname = splitext(bname)[0]
        if not bname: bname = "datablock"
        if len(params.image_range) == 1 and params.image_range[0] is not None:
          ext = "_{0}_{1}.json".format(*params.image_range[0])
        else:
          ext = "_sliced.json"
        output_datablocks_filename = bname + ext
      print('Saving sliced datablocks to {0}'.format(
        output_datablocks_filename))

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

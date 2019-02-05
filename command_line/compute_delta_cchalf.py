#!/usr/bin/env python
#
# dials.compute_delta_cchalf
#
#  Copyright (C) 2018 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#
from __future__ import absolute_import, division, print_function
import sys
from iotbx.reflection_file_reader import any_reflection_file
import matplotlib
matplotlib.use('Agg')
import dials.util
from matplotlib import pylab
from matplotlib import cm
from math import sqrt, floor
from cctbx import miller
from cctbx import crystal
from collections import defaultdict
from dials.algorithms.statistics.delta_cchalf import PerImageCChalfStatistics
from dials.array_family import flex
from dxtbx.model.experiment_list import ExperimentList
from libtbx.phil import parse
from dials.util import Sorry
import collections
import logging
from dials.util.exclude_images import exclude_image_ranges_for_scaling
from dials.util.multi_dataset_handling import select_datasets_on_ids
from dials.algorithms.scaling.scaling_library import set_image_ranges_in_scaling_models

logger = logging.getLogger('dials.command_line.compute_delta_cchalf')

help_message = '''

This program computes the delta cchalf excluding images

'''

# Set the phil scope
phil_scope = parse('''

  input {

    mtzfile = None
      .type = str
      .help = "We can also import an MTZ file"

  }

  mode = *dataset image_group
    .type = choice
    .help = "Perform analysis on whole datasets or batch groups"

  group_size = 10
    .type = int(value_min=1)
    .help = "The number of images to group together when calculating delta"
            "cchalf in image_group mode"

  output {

    experiments = "filtered_experiments.json"
      .type = str
      .help = "The filtered experiments file"

    reflections = "filtered_reflections.pickle"
      .type = str
      .help = "The filtered reflections file"

    table = "delta_cchalf.dat"
      .type = str
      .help = "A file with delta cchalf values"
  }

  nbins = 10
    .type = int(value_min=1)
    .help = "The number of resolution bins to use"

  dmin = None
    .type = float
    .help = "The maximum resolution"

  dmax = None
    .type = float
    .help = "The minimum resolution"

  stdcutoff = 4.0
    .type = float
    .help = "Datasets with a delta cc half below (mean - stdcutoff*std) are removed"

  output {

    log = 'dials.compute_delta_cchalf.log'
      .type = str
      .help = "The log filename"

    debug_log = 'dials.compute_delta_cchalf.debug.log'
      .type = str
      .help = "The debug log filename"
  }
include scope dials.util.multi_dataset_handling.phil_scope
''', process_includes=True)



class Script(object):
  '''A class for running the script.'''

  def __init__(self, params, experiments, reflections):
    '''Initialise the script.'''
    self.experiments = experiments
    self.reflections = reflections
    self.params = params
    # Set up a named tuple
    self.DataRecord = collections.namedtuple("DataRecord", (
      "unit_cell", "space_group", "miller_index", "dataset",
      "intensity", "variance", "identifiers", "images"))

  def prepare_data(self):
    if self.params.mode == 'image_group':
      for exp in self.experiments:
        if not exp.scan:
          raise Sorry("Cannot use mode=image_group with scanless experiments")

    if len(self.experiments) > 0 and len(self.reflections) == 1:
      data = self.read_experiments(self.experiments, self.reflections[0])
    elif len(self.experiments) > 0 and len(self.experiments) == len(self.reflections):
      # need to join together reflections
      joint_table = flex.reflection_table()
      for table in self.reflections:
        joint_table.extend(table)
      self.reflections = [joint_table]
      data = self.read_experiments(self.experiments, self.reflections[0])
    elif self.params.input.mtzfile is not None:
      data = self.read_mtzfile(self.params.input.mtzfile)
    else:
      self.parser.print_help()
      return SystemExit
    return data

  def run(self):
    '''Execute the script.'''
    from time import time

    data = self.prepare_data()

    # Create the statistics object
    statistics = PerImageCChalfStatistics(
      data.miller_index,
      data.identifiers,
      data.dataset,
      data.images,
      data.intensity,
      data.variance,
      data.unit_cell,
      data.space_group,
      nbins = self.params.nbins,
      dmin  = self.params.dmin,
      dmax  = self.params.dmax,
      mode = self.params.mode,
      image_group = self.params.group_size)

    if self.params.mode == 'image_group':
      self.image_group_to_id = statistics.image_group_to_id
      self.image_group_to_image_range = statistics.image_group_to_image_range

    # Print out the datasets in order of delta cc 1/2
    self.delta_cchalf_i = statistics.delta_cchalf_i()
    datasets = list(self.delta_cchalf_i.keys())
    sorted_index = sorted(range(len(datasets)), key=lambda x: self.delta_cchalf_i[datasets[x]])
    for i in sorted_index:
      logger.info("Dataset: %d, Delta CC 1/2: %.3f" % (datasets[i], 100*self.delta_cchalf_i[datasets[i]]))

    # Write a text file with delta cchalf values
    self.write_delta_cchalf_file(datasets, self.delta_cchalf_i, self.params)

    # Remove datasets based on delta cc1/2
    if len(self.experiments) > 0 and len(self.reflections) == 1:
      filtered_reflections = self.remove_datasets_below_cutoff(self.experiments,
        self.reflections[0], self.params, self.delta_cchalf_i)
      self.reflections = [filtered_reflections]



  def write_delta_cchalf_file(self, datasets, delta_cchalf_i, params):
    '''
    Write values to file

    '''
    logger.info("Writing table to %s" % params.output.table)
    with open(params.output.table, "w") as outfile:
      sorted_index = sorted(range(len(datasets)), key=lambda x: delta_cchalf_i[datasets[x]])
      for i in sorted_index:
        outfile.write("%d %f\n" % (datasets[i], 100*delta_cchalf_i[datasets[i]]))

  def read_experiments(self, experiments, reflections):
    '''
    Get information from experiments and reflections

    '''

    # Get space group and unit cell
    space_group = None
    unit_cell = []
    exp_identifiers = []
    for e in experiments:
      if space_group is None:
        space_group = e.crystal.get_space_group()
      else:
        assert space_group.type().number() == e.crystal.get_space_group().type().number()
      unit_cell.append(e.crystal.get_unit_cell())
      exp_identifiers.append(e.identifier)
    # get a list of the ids from the reflection table corresponding to exp_ids
    identifiers = []
    for expit in exp_identifiers:
      for k in reflections.experiment_identifiers().keys():
        if reflections.experiment_identifiers()[k] == expit:
          identifiers.append(k)
          break

    # Selection of reflections
    selection = ~(reflections.get_flags(reflections.flags.bad_for_scaling, all=False))
    outliers = reflections.get_flags(reflections.flags.outlier_in_scaling)
    reflections = reflections.select(selection & ~outliers)

    # Scale factor
    inv_scale_factor = reflections['inverse_scale_factor']
    selection = inv_scale_factor > 0
    reflections = reflections.select(selection)
    inv_scale_factor = reflections['inverse_scale_factor']

    # Get the reflection data
    index = reflections['id']

    miller_index = reflections['miller_index']
    intensity = reflections['intensity.scale.value'] / inv_scale_factor
    variance = reflections['intensity.scale.variance'] / inv_scale_factor**2
    # calculate image number of observation (e.g 0.0 <= z < 1.0), image = 1
    images = flex.floor(reflections['xyzobs.px.value'].parts()[2]).iround() + 1
    # Get the MTZ file
    return self.DataRecord(
      unit_cell    = unit_cell,
      space_group  = space_group,
      miller_index = miller_index,
      dataset      = index,
      intensity    = intensity,
      variance     = variance,
      identifiers  = identifiers,
      images = images)

  def read_mtzfile(self, filename):
    '''
    Read the mtz file

    '''
    # Read the mtz file
    reader = any_reflection_file(filename)

    # Get the columns as miller arrays
    miller_arrays = reader.as_miller_arrays(merge_equivalents=False)

    # Select the desired columns
    intensities = None
    batches = None
    for array in miller_arrays:
      if array.info().labels == ['I', 'SIGI']:
        intensities = array
      if array.info().labels == ['BATCH']:
        batches = array
    assert intensities is not None
    assert batches is not None
    assert len(batches.data()) == len(intensities.data())

    # Get the unit cell and space group
    unit_cell = intensities.unit_cell()
    space_group = intensities.crystal_symmetry().space_group()

    # The reflection data
    miller_index = intensities.indices()
    batch = batches.data()
    intensity = intensities.data()
    variance = intensities.sigmas()**2

    # Create unit cell list
    min_batch = min(batch)
    dataset = batch - min_batch
    num_datasets = max(dataset)+1
    unit_cell_list = [unit_cell for i in range(num_datasets)]

    # Get the MTZ file
    return self.DataRecord(
      unit_cell    = unit_cell_list,
      space_group  = space_group,
      miller_index = miller_index,
      dataset      = dataset,
      intensity    = intensity,
      variance     = variance)

  def remove_datasets_below_cutoff(self, experiments, reflections,
      params, delta_cchalf_i):
    '''
    Write the experiments and reflections

    '''
    X = list(delta_cchalf_i.keys())
    Y = list(delta_cchalf_i.values())
    mean = sum(Y) / len(Y)
    sdev = sqrt(sum((yy-mean)**2 for yy in Y)/len(Y))
    logger.info("\nmean delta_cc_half %s" % (mean*100))
    logger.info("stddev delta_cc_half %s" % (sdev*100))
    cutoff_value = mean - params.stdcutoff*sdev
    logger.info("cutoff value: %s \n" % (cutoff_value*100))
    datasets_to_remove = []
    if params.mode == 'dataset':
      datasets_to_remove = []
      for x in sorted(delta_cchalf_i.keys()):
        y = delta_cchalf_i[x]
        if y < cutoff_value:
          logger.info("Removing dataset %d" % x)
          datasets_to_remove.append(reflections.experiment_identifiers()[x])
      output_reflections = reflections.remove_on_experiment_identifiers(datasets_to_remove)
      experiments.remove_on_experiment_identifiers(datasets_to_remove)
      output_reflections.assert_experiment_identifiers_are_consistent(experiments)

    elif params.mode == 'image_group':
      datasets_to_remove = []
      exclude_images = []

      for x in sorted(delta_cchalf_i.keys()):
        y = delta_cchalf_i[x]
        if y < cutoff_value:
          exp_id = self.image_group_to_id[x] #numerical id
          identifier = reflections.experiment_identifiers()[exp_id]
          image_range = self.image_group_to_image_range[x]
          datasets_to_remove.append(exp_id)
          logger.info("Removing image range %s from experiment %s",
            image_range, identifier)
          exclude_images.append(
            [identifier+':'+str(image_range[0])+':'+str(image_range[1])])

      # Now remove individual batches
      if -1 in reflections['id']:
        reflections = reflections.select(reflections['id'] != -1)
      reflection_list = reflections.split_by_experiment_id()
      reflection_list, experiments = exclude_image_ranges_for_scaling(
        reflection_list, experiments, exclude_images)
      #if a whole experiment has been excluded: need to remove it here
      experiments_to_delete = []
      for exp in experiments:
        if not exp.scan.get_valid_image_ranges(exp.identifier): #if all removed above
          experiments_to_delete.append(exp.identifier)
      if experiments_to_delete:
        experiments, reflection_list = select_datasets_on_ids(
          experiments, reflection_list, exclude_datasets=experiments_to_delete)
      assert len(reflection_list) == len(experiments)
      experiments = set_image_ranges_in_scaling_models(experiments)

      output_reflections = flex.reflection_table()
      for r in reflection_list:
        output_reflections.extend(r)
    return output_reflections

    # Write the experiments and reflections to file
  def write_experiments_and_reflections(self):
    if len(self.experiments) > 0 and len(self.reflections) == 1:
      self.write_reflections(self.reflections[0], self.params.output.reflections)
      self.write_experiments(self.experiments, self.params.output.experiments)

  def write_reflections(self, reflections, filename):
    ''' Save the reflections to file. '''
    logger.info('Saving %d reflections to %s' % (len(reflections), filename))
    reflections.as_pickle(filename)

  def write_experiments(self, experiments, filename):
    ''' Save the profile model parameters. '''
    from dxtbx.model.experiment_list import ExperimentListDumper
    logger.info('Saving the experiments to %s' % filename)
    dump = ExperimentListDumper(experiments)
    with open(filename, "w") as outfile:
      outfile.write(dump.as_json())

  def plot_data(self):
    # Make a plot of delta cc 1/2
    from matplotlib import pylab
    fig, ax = pylab.subplots()
    ax.hist(self.delta_cchalf_i.values())
    ax.set_xlabel("Delta CC 1/2")
    fig.savefig("plot1.png")
    #pylab.show()

    X = list(self.delta_cchalf_i.keys())
    Y = list(self.delta_cchalf_i.values())
    fig, ax = pylab.subplots()
    ax.plot(X, Y)
    ax.set_xlabel("Dataset number")
    ax.set_ylabel("Delta CC 1/2")
    fig.savefig("plot2.png")
  #pylab.show()

def run(args=None, phil=phil_scope):
  import dials.util.log
  from dials.util.options import OptionParser
  from dials.util.options import flatten_reflections
  from dials.util.options import flatten_experiments

  usage = "dials.compute_delta_cchalf [options] scaled_experiments.json scaled.pickle"

  parser = OptionParser(
    usage=usage,
    phil=phil,
    epilog=__doc__,
    read_experiments=True,
    read_reflections=True,
    check_format=False)

  params, _ = parser.parse_args(args=args, show_diff_phil=False)

  dials.util.log.config(info=params.output.log, debug=params.output.debug_log)

  expermients = flatten_experiments(params.input.experiments)
  reflections = flatten_reflections(params.input.reflections)

  script = Script(params, expermients, reflections)
  script.run()
  script.write_experiments_and_reflections()
  script.plot_data()

if __name__ == '__main__':
  with dials.util.show_mail_on_error():
    run()

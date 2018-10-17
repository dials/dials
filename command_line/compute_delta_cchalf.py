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
import collections
import logging

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

  output {

    experiments = "filtered_experiments.json"
      .type = str
      .help = "The filtered experiments file"

    reflections = "filtered_reflections.pickle"
      .type = str
      .help = "The filtered reflections file"
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

''', process_includes=True)


class Script(object):
  '''A class for running the script.'''

  def __init__(self):
    '''Initialise the script.'''
    from dials.util.options import OptionParser
    import libtbx.load_env

    # The script usage
    usage = "usage: %s [options] [param.phil] " % libtbx.env.dispatcher_name

    # Initialise the base class
    self.parser = OptionParser(
      usage=usage,
      phil=phil_scope,
      epilog=help_message,
      read_experiments=True,
      read_reflections=True,
      check_format=False)

  def run(self):
    '''Execute the script.'''
    from dials.array_family import flex
    from dials.util.options import flatten_experiments
    from dials.util.options import flatten_reflections
    from time import time
    from dials.util import log
    from libtbx.utils import Sorry

    # Parse the command line
    params, options = self.parser.parse_args(show_diff_phil=False)

    # Configure the logging
    log.config(
      info=params.output.log,
      debug=params.output.debug_log)

    from dials.util.version import dials_version
    logger.info(dials_version())

    # Log the diff phil
    diff_phil = self.parser.diff_phil.as_str()
    if diff_phil is not '':
      logger.info('The following parameters have been modified:\n')
      logger.info(diff_phil)

    # Setup a named tuple
    self.DataRecord = collections.namedtuple("DataRecord", (
      "unit_cell",
      "space_group",
      "miller_index",
      "dataset",
      "intensity",
      "variance"))

    # Ensure we have an experiment list
    experiments = flatten_experiments(params.input.experiments)
    reflections = flatten_reflections(params.input.reflections)
    if len(experiments) > 0:
      assert len(reflections) == 1
      data = self.read_experiments(experiments, reflections[0])
    elif params.input.mtzfile is not None:
      data = self.read_mtzfile(params.input.mtzfile)
    else:
      self.parser.print_help()
      return

    # Create the statistics object
    statistics = PerImageCChalfStatistics(
      data.miller_index,
      data.dataset,
      data.intensity,
      data.variance,
      data.unit_cell,
      data.space_group,
      nbins = params.nbins,
      dmin  = params.dmin,
      dmax  = params.dmax)

    # Print out the datasets in order of delta cc 1/2
    delta_cchalf_i = statistics.delta_cchalf_i()
    datasets = list(delta_cchalf_i.keys())
    sorted_index = sorted(range(len(datasets)), key=lambda x: delta_cchalf_i[datasets[x]])
    for i in sorted_index:
      print("Dataset: %d, Delta CC 1/2: %.3f" % (datasets[i], 100*delta_cchalf_i[datasets[i]]))

    # Remove datasets based on delta cc1/2
    if len(experiments) > 0:
      self.write_experiments_and_reflections(
        experiments,
        reflections[0],
        params,
        delta_cchalf_i)

    # Make a plot of delta cc 1/2
    from matplotlib import pylab
    fig, ax = pylab.subplots()
    ax.hist(delta_cchalf_i.values())
    ax.set_xlabel("Delta CC 1/2")
    pylab.show()

    X = list(delta_cchalf_i.keys())
    Y = list(delta_cchalf_i.values())
    fig, ax = pylab.subplots()
    ax.plot(X, Y)
    ax.set_xlabel("Dataset number")
    ax.set_ylabel("Delta CC 1/2")
    pylab.show()

  def read_experiments(self, experiments, reflections):
    '''
    Get information from experiments and reflections

    '''

    # Get space group and unit cell
    space_group = None
    unit_cell = []
    for e in experiments:
      if space_group is None:
        space_group = e.crystal.get_space_group()
      else:
        assert space_group.type().number() == e.crystal.get_space_group().type().number()
      unit_cell.append(e.crystal.get_unit_cell())

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

    # Get the MTZ file
    return self.DataRecord(
      unit_cell    = unit_cell,
      space_group  = space_group,
      miller_index = miller_index,
      dataset      = index,
      intensity    = intensity,
      variance     = variance)

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

  def write_experiments_and_reflections(self,
                                        experiments,
                                        reflections,
                                        params,
                                        delta_cchalf_i):
    '''
    Write the experiments and reflections

    '''
    X = list(delta_cchalf_i.keys())
    Y = list(delta_cchalf_i.values())
    mean = sum(Y) / len(Y)
    sdev = sqrt(sum((yy-mean)**2 for yy in Y)/len(Y))
    output_experiments = ExperimentList()
    output_reflections = reflections
    print("\nmean delta_cc_half %s" % (mean*100))
    print("stddev delta_cc_half %s" % (sdev*100))
    cutoff_value = mean - params.stdcutoff*sdev
    print("cutoff value: %s \n" % (cutoff_value*100))
    for x in sorted(delta_cchalf_i.keys()):
      y = delta_cchalf_i[x]
      if y < cutoff_value:
        print("Removing dataset %d" % x)
        output_reflections.del_selected(output_reflections['id'] == x)
        selection = output_reflections['id'] > x
        indices = output_reflections['id'].select(selection)
        indices -= 1
        output_reflections['id'].set_selected(selection, indices)
      else:
        output_experiments.append(experiments[x])
    assert max(output_reflections['id']) < len(output_experiments)

    # Write the experiments and reflections to file
    self.write_reflections(output_reflections, params.output.reflections)
    self.write_experiments(output_experiments, params.output.experiments)

  def write_reflections(self, reflections, filename):
    ''' Save the reflections to file. '''
    print('Saving %d reflections to %s' % (len(reflections), filename))
    reflections.as_pickle(filename)

  def write_experiments(self, experiments, filename):
    ''' Save the profile model parameters. '''
    from dxtbx.model.experiment_list import ExperimentListDumper
    print('Saving the experiments to %s' % filename)
    dump = ExperimentListDumper(experiments)
    with open(filename, "w") as outfile:
      outfile.write(dump.as_json())

if __name__ == '__main__':
  from dials.util import halraiser
  try:
    script = Script()
    script.run()
  except Exception as e:
    halraiser(e)

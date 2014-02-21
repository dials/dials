#!/usr/bin/env python
#
# dials.process.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division
from dials.util.script import ScriptRunner


class Script(ScriptRunner):
  '''A class for running the script.'''

  def __init__(self):
    '''Initialise the script.'''

    # The script usage
    usage = "usage: %prog [options] [param.phil] datablock.json"

    # Initialise the base class
    ScriptRunner.__init__(self, usage=usage)

    # Add a verbosity option
    self.config().add_option(
      "-v",
      dest="verbosity",
      action="count", default=1,
      help="set verbosity level; -vv gives verbosity level 2")

    # Add an options for spot finder output
    self.config().add_option(
      "--strong-filename",
      dest="strong_filename",
      type="str", default=None,
      help="The output filename for found spots")

  def main(self, params, options, args):
    '''Execute the script.'''
    from dials.util.command_line import Importer
    from dials.util.command_line import Command
    from time import time

    # Save the options
    self.options = options
    self.params = params

    st = time()

    # Preamble stuff
    print '*' * 80
    print ''
    print '                       mmmm   mmmmm    mm   m       mmmm            '
    print '                       #   "m   #      ##   #      #"   "           '
    print '                      m#mm  #   #     #  #  #      "#mmm            '
    print '                       #    #   #     #mm#  #          "#           '
    print '                       #mmm"  mm#mm  #    # #mmmmm "mmm#"           '
    print ''
    print 'Launching dials.process'
    print ''
    print 'The following tasks will be performed:'
    print ' 1) Strong spots will be found (dials.find_spots)'
    print ' 2) The strong spots will be indexed (dials.index)'
    print ' 3) A profile model will be created (dials.create_profile_model)'
    print ' 4) The reflections will be integrated (dials.integrate)'
    print ''
    print 'Please be patient, this may take a few minutes'
    print ''
    print '*' * 80
    print ''

    # Import stuff
    Command.start('Importing datablocks')
    importer = Importer(args, include=["images", "datablocks"])
    assert(len(importer.datablocks) == 1)
    datablock = importer.datablocks[0]
    Command.end('Imported datablocks')

    # Check the unhandled arguments
    if len(importer.unhandled_arguments) > 0:
      print '-' * 80
      print 'The following command line arguments weren\'t handled'
      for arg in importer.unhandled_arguments:
        print '  ' + arg
    print ''

    # Do the processing
    observed = self.find_spots(datablock)
    experiments, indexed = self.index(datablock, observed)
    profile = self.create_profile_model(experiments, indexed)
    integrated = self.integrate(experiments, profile, indexed)

    # Total Time
    print ""
    print "Total Time Taken = %f seconds" % (time() - st)

  def find_spots(self, datablock):
    from time import time
    from dials.array_family import flex
    from dials.util.command_line import Command
    st = time()

    print '*' * 80
    print 'Finding Strong Spots'
    print '*' * 80

    # Find the strong spots
    observed = flex.reflection_table.from_observations(datablock)

    # Save the reflections to file
    print '\n' + '-' * 80
    if self.options.strong_filename:
      Command.start('Saving {0} reflections to {1}'.format(
          len(observed), self.options.strong_filename))
      observed.as_pickle(self.options.strong_filename)
      Command.end('Saved {0} observed to {1}'.format(
          len(observed), self.options.strong_filename))

    print ''
    print 'Time Taken = %f seconds' % (time() - st)
    return observed

  def index(self, datablock, reflections):
    from dials.algorithms.indexing.indexer2 import master_phil_scope
    from libtbx.phil import command_line
    from time import time
    st = time()

    print '*' * 80
    print 'Indexing Strong Spots'
    print '*' * 80

    imagesets = datablock.extract_imagesets()
    if len(imagesets) > 1:
      raise RuntimeError("Only one imageset can be processed at a time")
    imageset = imagesets[0]

    cmd_line = command_line.argument_interpreter(master_params=master_phil_scope)
    working_phil = cmd_line.process_and_fetch(args=[])
    working_phil.show()

    gonio = imageset.get_goniometer()
    detector = imageset.get_detector()
    scan = imageset.get_scan()
    beam = imageset.get_beam()
    print detector
    print scan
    print gonio
    print beam

    params = working_phil.extract()
    params.refinement.reflections.use_all_reflections=True
    if params.method == "fft3d":
      from dials.algorithms.indexing.fft3d import indexer_fft3d as indexer
    elif params.method == "fft1d":
      from dials.algorithms.indexing.fft1d import indexer_fft1d as indexer
    elif params.method == "real_space_grid_search":
      from dials.algorithms.indexing.real_space_grid_search \
           import indexer_real_space_grid_search as indexer
    idxr = indexer(reflections, imageset, params=params)
    idxr.index()

    from dials.array_family import flex
    from dials.model.experiment.experiment_list import ExperimentListFactory
    indexed = flex.reflection_table.from_pickle("indexed.pickle")
    experiments = ExperimentListFactory.from_json_file("experiments.json")

    print ''
    print 'Time Taken = %f seconds' % (time() - st)
    return experiments, indexed

  def create_profile_model(self, experiments, reflections):
    from dials.algorithms.profile_model.profile_model import ProfileModel
    from time import time
    from dials.util.command_line import Command
    from math import pi
    st = time()

    print '*' * 80
    print 'Creating Profile Model'
    print '*' * 80

    print 'Starting with %d reflections' % len(reflections)

    from dials.array_family import flex
    Command.start('Removing invalid coordinates')
    xyz = reflections['xyzcal.mm']
    mask = flex.bool([x == (0, 0, 0) for x in xyz])
    reflections.del_selected(mask)
    Command.end('Removed invalid coordinates, %d remaining' % len(reflections))

    # Create the profile model
    profile_model = ProfileModel(experiments[0], reflections)
    sigma_b = profile_model.sigma_b() * 180.0 / pi
    sigma_m = profile_model.sigma_m() * 180.0 / pi
    print 'Sigma B: %f' % sigma_b
    print 'Sigma M: %f' % sigma_m

    # Write the parameters
    from dials.framework.registry import Registry
    registry = Registry()

    # Get the parameters
    params = registry.config().params()
    params.shoebox.sigma_b = sigma_b
    params.shoebox.sigma_m = sigma_m

    print ''
    print 'Time Taken = %f seconds' % (time() - st)
    return (sigma_b, sigma_m)

  def integrate(self, experiments, profile, indexed):
    from dials.algorithms.integration import Integrator
    from time import time
    from dials.framework.registry import Registry
    registry = Registry()
    assert(registry.params().shoebox.sigma_b > 0)
    assert(registry.params().shoebox.sigma_m > 0)

    st = time()

    print '*' * 80
    print 'Integrating Reflections'
    print '*' * 80

    # Get the integrator from the input parameters
    print 'Configurating integrator from input parameters'
    integrate = Integrator(self.params.shoebox.n_sigma,
                           self.params.shoebox.n_blocks,
                           self.params.integration.filter.by_zeta)

    # Intregate the sweep's reflections
    print 'Integrating reflections'
    reflections = integrate(experiments, reference=indexed, extracted=None)

    print ''
    print 'Time Taken = %f seconds' % (time() - st)
    return reflections


if __name__ == '__main__':
  script = Script()
  script.run()

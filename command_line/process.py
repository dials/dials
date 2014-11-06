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

help_message = '''

This is the main dials program. Given a set of images, the data will be
processed and a list of integrated reflections will be given as output. More
specialised help can be seen by looking at the command line programs that
perform the individual processing steps. Looks at the following for more
details:

  dials.import
  dials.find_spots
  dials.index
  dials.refine
  dials.integrate
  dials.export_mtz

This program will do the following:

  First the image data will be imported into a datablock. Strong spots will then
  be found on all the images of the datablock. These strong spots will then be
  indexed (the indexing step also includes some static centroid refinement). The
  experimental geometry will then be refined and the reflections integrated.
  Finally, the integrated reflections will be exported to an MTZ file which can
  be input into Aimless to be scaled.

  Examples:

    dials.process images*.cbf

    dials.process datablock.json

    find . -name "images*.cbf" | dials.process

'''

class Script(object):
  '''A class for running the script.'''

  def __init__(self):
    '''Initialise the script.'''
    from dials.util.options import OptionParser
    from libtbx.phil import parse
    import libtbx.load_env

    phil_scope = parse('''
      output {
        datablock_filename = datablock.json
          .type = str
          .help = "The filename for output datablock"

        strong_filename = strong.pickle
          .type = str
          .help = "The filename for strong reflections from spot finder output."

        indexed_filename = indexed.pickle
          .type = str
          .help = "The filename for indexed reflections."

        refined_experiments_filename = refined_experiments.json
          .type = str
          .help = "The filename for saving refined experimental models"

        integrated_filename = integrated.pickle
          .type = str
          .help = "The filename for final integrated reflections."

        profile_filename = profile.phil
          .type = str
          .help = "The filename for output reflection profile parameters"

        mtz_filename = integrated.mtz
          .type = str
          .help = "The filename for output mtz"
      }

      include scope dials.algorithms.peak_finding.spotfinder_factory.phil_scope
      include scope dials.algorithms.indexing.indexer.index_only_phil_scope
      include scope dials.algorithms.refinement.refiner.phil_scope
      include scope dials.algorithms.integration.integrator.phil_scope
      include scope dials.algorithms.profile_model.factory.phil_scope
      include scope dials.algorithms.spot_prediction.reflection_predictor.phil_scope

    ''', process_includes=True)

    # The script usage
    usage = "usage: %s [options] [param.phil] datablock.json" % libtbx.env.dispatcher_name

    # Create the parser
    self.parser = OptionParser(
      usage=usage,
      phil=phil_scope,
      epilog=help_message,
      read_datablocks=True,
      read_datablocks_from_images=True)

  def run(self):
    '''Execute the script.'''
    from dials.util.command_line import Command
    from dials.util.options import flatten_datablocks
    from time import time
    import sys

    # Parse the command line
    params, options = self.parser.parse_args(show_diff_phil=True)
    datablocks = flatten_datablocks(params.input.datablock)

    # Save the options
    self.options = options
    self.params = params

    st = time()

    # Import stuff
    if len(datablocks) == 0:
      self.parser.print_help()
      return
    elif len(datablocks) > 1:
      raise RuntimeError('Only 1 datablock can be processed at a time.')
    datablock = datablocks[0]

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
    print ' 3) The model will be further refined (dials.refine)'
    print ' 4) The reflections will be integrated (dials.integrate)'
    print ' 5) The data will be exported as MTZ (dials.export_mtz)'
    print ''
    print 'Please be patient, this may take a few minutes'
    print ''
    print '*' * 80
    print ''
    print 'Command-line: %s' % (' '.join(sys.argv[1:]))[:65]
    print ''
    print '*' * 80
    print ''

    if self.params.output.datablock_filename:
      from dxtbx.datablock import DataBlockDumper
      dump = DataBlockDumper(datablock)
      dump.as_json(self.params.output.datablock_filename)

    # Do the processing
    observed = self.find_spots(datablock)
    experiments, indexed = self.index(datablock, observed)
    experiments = self.refine(experiments, indexed)
    integrated = self.integrate(experiments, indexed)
    mtz = self.mtz(integrated, experiments)
    mtz.show_summary()

    # Total Time
    print ""
    print "Total Time Taken = %f seconds" % (time() - st)

  def find_spots(self, datablock):
    from time import time
    from dials.array_family import flex
    st = time()

    print '*' * 80
    print 'Finding Strong Spots'
    print '*' * 80

    # Find the strong spots
    observed = flex.reflection_table.from_observations(datablock, self.params)

    # Save the reflections to file
    print '\n' + '-' * 80
    if self.params.output.strong_filename:
      from dials.util.command_line import Command
      Command.start('Saving {0} reflections to {1}'.format(
          len(observed), self.params.output.strong_filename))
      observed.as_pickle(self.params.output.strong_filename)
      Command.end('Saved {0} observed to {1}'.format(
          len(observed), self.params.output.strong_filename))

    print ''
    print 'Time Taken = %f seconds' % (time() - st)
    return observed

  def index(self, datablock, reflections):
    from time import time
    import copy
    st = time()

    print '*' * 80
    print 'Indexing Strong Spots'
    print '*' * 80

    imagesets = datablock.extract_imagesets()

    params = copy.deepcopy(self.params)
    # don't do scan-varying refinement during indexing
    params.refinement.parameterisation.crystal.scan_varying = False
    if params.indexing.method == "fft3d":
      from dials.algorithms.indexing.fft3d import indexer_fft3d as indexer
    elif params.indexing.method == "fft1d":
      from dials.algorithms.indexing.fft1d import indexer_fft1d as indexer
    elif params.indexing.method == "real_space_grid_search":
      from dials.algorithms.indexing.real_space_grid_search \
           import indexer_real_space_grid_search as indexer
    idxr = indexer(reflections, imagesets, params=params)

    indexed = idxr.refined_reflections
    experiments = idxr.refined_experiments

    if self.params.output.indexed_filename:
      from dials.util.command_line import Command
      Command.start('Saving {0} reflections to {1}'.format(
          len(indexed), self.params.output.indexed_filename))
      indexed.as_pickle(self.params.output.indexed_filename)
      Command.end('Saved {0} reflections to {1}'.format(
          len(indexed), self.params.output.indexed_filename))

    print ''
    print 'Time Taken = %f seconds' % (time() - st)
    return experiments, indexed

  def refine(self, experiments, centroids):
    from dials.algorithms.refinement import RefinerFactory
    from time import time
    st = time()

    print '*' * 80
    print 'Refining Model'
    print '*' * 80

    refiner = RefinerFactory.from_parameters_data_experiments(
      self.params, centroids, experiments)

    refiner.run()
    experiments = refiner.get_experiments()

    # Dump experiments to disk
    if self.params.output.refined_experiments_filename:
      from dxtbx.model.experiment.experiment_list import ExperimentListDumper
      dump = ExperimentListDumper(experiments)
      dump.as_json(self.params.output.refined_experiments_filename)

    print ''
    print 'Time Taken = %f seconds' % (time() - st)

    return experiments

  def integrate(self, experiments, indexed):
    from time import time
    from dials.util.command_line import Command

    st = time()

    print '*' * 80
    print 'Integrating Reflections'
    print '*' * 80


    from dials.array_family import flex
    Command.start('Removing invalid coordinates')
    xyz = indexed['xyzcal.mm']
    mask = flex.bool([x == (0, 0, 0) for x in xyz])
    indexed.del_selected(mask)
    Command.end('Removed invalid coordinates, %d remaining' % len(indexed))

    # Get the integrator from the input parameters
    print 'Configurating integrator from input parameters'
    from dials.algorithms.profile_model.factory import ProfileModelFactory
    from dials.algorithms.integration.integrator import IntegratorFactory
    from dials.array_family import flex

    # Compute the profile model
    # Predict the reflections
    # Match the predictions with the reference
    # Create the integrator
    profile_model = ProfileModelFactory.create(self.params, experiments, indexed)
    print ""
    print "=" * 80
    print ""
    print "Predicting reflections"
    print ""
    predicted = profile_model.predict_reflections(
      experiments,
      dmin=self.params.prediction.dmin,
      dmax=self.params.prediction.dmax,
      margin=self.params.prediction.margin,
      force_static=self.params.prediction.force_static)
    predicted.match_with_reference(indexed)
    print ""
    integrator = IntegratorFactory.create(self.params, experiments, profile_model, predicted)

    # Integrate the reflections
    reflections = integrator.integrate()

    if self.params.output.integrated_filename:
      # Save the reflections
      Command.start('Saving {0} reflections to {1}'.format(
          len(reflections), self.params.output.integrated_filename))
      reflections.as_pickle(self.params.output.integrated_filename)
      Command.end('Saved {0} reflections to {1}'.format(
          len(reflections), self.params.output.integrated_filename))

    if self.params.output.profile_filename:
      with open(self.params.output.profile_filename, "w") as outfile:
        outfile.write(profile_model.dump().as_str())

    print ''
    print 'Time Taken = %f seconds' % (time() - st)
    return reflections

  def mtz(self, integrated, experiments):
    from dials.util.export_mtz import export_mtz
    from time import time
    st = time()
    print '*' * 80
    print 'Exporting measurements to', self.params.output.mtz_filename
    print '*' * 80
    m = export_mtz(integrated, experiments, self.params.output.mtz_filename)
    print ''
    print 'Time Taken = %f seconds' % (time() - st)
    return m

if __name__ == '__main__':
  from dials.util import halraiser
  try:
    script = Script()
    script.run()
  except Exception as e:
    halraiser(e)

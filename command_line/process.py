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

    find . -name "images*.cbf" | dials.process -i

'''

class Script(object):
  '''A class for running the script.'''

  def __init__(self):
    '''Initialise the script.'''
    from dials.util.options import OptionParser
    from libtbx.phil import parse

    phil_scope = parse('''

      strong_filename = None
        .type = str
        .help = "The filename for strong reflections from spot finder output."

      include scope dials.algorithms.peak_finding.spotfinder_factory.phil_scope
      indexing {
        include scope dials.algorithms.indexing.indexer.master_phil_scope
      }
      include scope dials.algorithms.refinement.refiner.phil_scope
      include scope dials.algorithms.integration.interface.phil_scope
      include scope dials.algorithms.profile_model.factory.phil_scope
      include scope dials.algorithms.spot_prediction.reflection_predictor.phil_scope

    ''', process_includes=True)

    # The script usage
    usage = "usage: %prog [options] [param.phil] datablock.json"

    # Create the parser
    self.parser = OptionParser(
      usage=usage,
      stdin_options=True,
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

    # Parse the command line
    params, options = self.parser.parse_args(show_diff_phil=True)
    datablocks = flatten_datablocks(params.input.datablock)

    # Save the options
    self.options = options
    self.params = params

    st = time()

    # Import stuff
    assert(len(datablocks) == 1)
    datablock = datablocks[0]

    from dxtbx.datablock import DataBlockDumper
    dump = DataBlockDumper(datablock)
    dump.as_json("datablock.json")

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
    from dials.util.command_line import Command
    st = time()

    print '*' * 80
    print 'Finding Strong Spots'
    print '*' * 80

    # Find the strong spots
    observed = flex.reflection_table.from_observations(datablock)

    # Save the reflections to file
    print '\n' + '-' * 80
    if self.params.strong_filename:
      Command.start('Saving {0} reflections to {1}'.format(
          len(observed), self.params.strong_filename))
      observed.as_pickle(self.params.strong_filename)
      Command.end('Saved {0} observed to {1}'.format(
          len(observed), self.params.strong_filename))

    observed.as_pickle("strong.pickle")

    print ''
    print 'Time Taken = %f seconds' % (time() - st)
    return observed

  def index(self, datablock, reflections):
    from time import time
    st = time()

    print '*' * 80
    print 'Indexing Strong Spots'
    print '*' * 80

    imagesets = datablock.extract_imagesets()

    #params = working_phil.extract()
    params = self.params.indexing
    if params.method == "fft3d":
      from dials.algorithms.indexing.fft3d import indexer_fft3d as indexer
    elif params.method == "fft1d":
      from dials.algorithms.indexing.fft1d import indexer_fft1d as indexer
    elif params.method == "real_space_grid_search":
      from dials.algorithms.indexing.real_space_grid_search \
           import indexer_real_space_grid_search as indexer
    idxr = indexer(reflections, imagesets, params=params)

    indexed = idxr.refined_reflections
    experiments = idxr.refined_experiments

    indexed.as_pickle("indexed.pickle")

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
    from dxtbx.model.experiment.experiment_list import ExperimentListDumper
    dump = ExperimentListDumper(experiments)
    dump.as_json("refined_experiments.json")

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
    from dials.algorithms.integration.interface import IntegratorFactory
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
    predicted = flex.reflection_table.from_predictions_multi(
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

    # Save the reflections
    reflections.as_pickle("integrated.pickle")
    with open("profile.phil", "w") as outfile:
      outfile.write(profile_model.dump().as_str())

    print ''
    print 'Time Taken = %f seconds' % (time() - st)
    return reflections

  def mtz(self, integrated, experiments):
    from dials.util.export_mtz import export_mtz
    from time import time
    st = time()
    print '*' * 80
    print 'Exporting measurements to integrated.mtz'
    print '*' * 80
    m = export_mtz(integrated, experiments, 'integrated.mtz')
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

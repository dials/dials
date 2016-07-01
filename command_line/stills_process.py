#!/usr/bin/env python
#
# LIBTBX_SET_DISPATCHER_NAME dials.stills_process

from __future__ import division
from libtbx.utils import Sorry

help_message = '''
DIALS script for processing still images. Import, index, refine, and integrate are all done for each image
seperately.
'''

from libtbx.phil import parse
control_phil_str = '''
  verbosity = 1
    .type = int(value_min=0)
    .help = "The verbosity level"

  output {
    output_dir = .
      .type = str
      .help = Directory output files will be placed
    logging_dir = None
      .type = str
      .help = Directory output log files will be placed
    datablock_filename = %s_datablock.json
      .type = str
      .help = The filename for output datablock
    strong_filename = %s_strong.pickle
      .type = str
      .help = The filename for strong reflections from spot finder output.
    indexed_filename = %s_indexed.pickle
      .type = str
      .help = The filename for indexed reflections.
    refined_experiments_filename = %s_refined_experiments.json
      .type = str
      .help = The filename for saving refined experimental models
    integrated_filename = %s_integrated.pickle
      .type = str
      .help = The filename for final integrated reflections.
    profile_filename = None
      .type = str
      .help = The filename for output reflection profile parameters
    integration_pickle = int-%d-%s.pickle
      .type = str
      .help = Filename for cctbx.xfel-style integration pickle files
  }

  mp {
    method = *multiprocessing sge lsf pbs
      .type = choice
      .help = "The multiprocessing method to use"
    nproc = 1
      .type = int(value_min=1)
      .help = "The number of processes to use."
  }
'''

dials_phil_str = '''

  include scope dials.algorithms.spot_finding.factory.phil_scope
  include scope dials.algorithms.indexing.indexer.index_only_phil_scope
  include scope dials.algorithms.refinement.refiner.phil_scope
  include scope dials.algorithms.integration.integrator.phil_scope
  include scope dials.algorithms.profile_model.factory.phil_scope
  include scope dials.algorithms.spot_prediction.reflection_predictor.phil_scope

  integration {
    summation {
      detector_gain = 1
        .type = float
        .help = Multiplier for variances after integration. See Leslie 1999.
    }
  }
'''

phil_scope = parse(control_phil_str + dials_phil_str, process_includes=True)

class Script(object):
  '''A class for running the script.'''

  def __init__(self):
    '''Initialise the script.'''
    from dials.util.options import OptionParser
    import libtbx.load_env

    # The script usage
    usage = "usage: %s [options] [param.phil] filenames" % libtbx.env.dispatcher_name

    # Create the parser
    self.parser = OptionParser(
      usage=usage,
      phil=phil_scope,
      epilog=help_message
      )

  def run(self):
    '''Execute the script.'''
    from dials.util import log
    from logging import info
    from time import time
    from libtbx.utils import Abort
    from libtbx import easy_mp
    import os, copy
    from dxtbx.datablock import DataBlockFactory

    # Parse the command line
    params, options, all_paths = self.parser.parse_args(show_diff_phil=False, return_unhandled=True)

    # Check we have some filenames
    if not all_paths:
      self.parser.print_help()
      return

    # Save the options
    self.options = options
    self.params = params

    st = time()

    # Configure logging
    log.config(
      params.verbosity,
      info='dials.process.log',
      debug='dials.process.debug.log')

    # Log the diff phil
    diff_phil = self.parser.diff_phil.as_str()
    if diff_phil is not '':
      info('The following parameters have been modified:\n')
      info(diff_phil)

    # Import stuff
    info("Loading files...")
    def do_import(filename):
      info("Loading %s"%os.path.basename(filename))
      try:
        datablocks = DataBlockFactory.from_json_file(filename)
      except ValueError:
        datablocks = DataBlockFactory.from_filenames([filename])
      if len(datablocks) == 0:
        raise Abort("Could not load %s"%filename)
      if len(datablocks) > 1:
        raise Abort("Got multiple datablocks from file %s"%filename)
      return datablocks[0]

    if len(all_paths) == 1:
      datablocks = [do_import(all_paths[0])]
    else:
      datablocks = easy_mp.parallel_map(
        func=do_import,
        iterable=all_paths,
        processes=params.mp.nproc,
        method=params.mp.method,
        preserve_order=True,
        preserve_exception_message=True)

    if len(datablocks) == 0:
      raise Abort('No datablocks specified')

    # Handle still imagesets by breaking them apart into multiple datablocks
    # Further handle single file still imagesets (like HDF5) by tagging each
    # frame using its index
    indices = []
    basenames = []
    split_datablocks = []
    for datablock in datablocks:
      for imageset in datablock.extract_imagesets():
        for i in xrange(len(imageset)):
          subset = imageset[i:i+1]
          split_datablocks.append(DataBlockFactory.from_imageset(subset)[0])
          indices.append(i)
          basenames.append(os.path.splitext(os.path.basename(subset.paths()[0]))[0])
    tags = []
    for i, basename in zip(indices, basenames):
      if basenames.count(basename) > 1:
        tags.append("%s_%d"%(basename, i))
      else:
        tags.append(basename)

    # Wrapper function
    def do_work(item):
      Processor(copy.deepcopy(params)).process_datablock(item[0], item[1])

    # Process the data
    easy_mp.parallel_map(
      func=do_work,
      iterable=zip(tags, split_datablocks),
      processes=params.mp.nproc,
      method=params.mp.method,
      preserve_order=True,
      preserve_exception_message=True)

     # Total Time
    info("")
    info("Total Time Taken = %f seconds" % (time() - st))

class Processor(object):
  def __init__(self, params):
    self.params = params

  def process_datablock(self, tag, datablock):
    import os
    s = tag

    # before processing, set output paths according to the templates
    if self.params.output.datablock_filename is not None and "%s" in self.params.output.datablock_filename:
      self.params.output.datablock_filename = os.path.join(self.params.output.output_dir, self.params.output.datablock_filename%("idx-" + s))
    if self.params.output.strong_filename is not None and "%s" in self.params.output.strong_filename:
      self.params.output.strong_filename = os.path.join(self.params.output.output_dir, self.params.output.strong_filename%("idx-" + s))
    if self.params.output.indexed_filename is not None and "%s" in self.params.output.indexed_filename:
      self.params.output.indexed_filename = os.path.join(self.params.output.output_dir, self.params.output.indexed_filename%("idx-" + s))
    if "%s" in self.params.output.refined_experiments_filename:
      self.params.output.refined_experiments_filename = os.path.join(self.params.output.output_dir, self.params.output.refined_experiments_filename%("idx-" + s))
    if "%s" in self.params.output.integrated_filename:
      self.params.output.integrated_filename = os.path.join(self.params.output.output_dir, self.params.output.integrated_filename%("idx-" + s))

    if self.params.output.datablock_filename:
      from dxtbx.datablock import DataBlockDumper
      dump = DataBlockDumper(datablock)
      dump.as_json(self.params.output.datablock_filename)

    # Do the processing
    try:
      observed = self.find_spots(datablock)
    except Exception, e:
      print "Error spotfinding", tag, str(e)
      return
    try:
      experiments, indexed = self.index(datablock, observed)
    except Exception, e:
      print "Couldn't index", tag, str(e)
      return
    try:
      experiments = self.refine(experiments, indexed)
    except Exception, e:
      print "Error refining", tag, str(e)
      return
    try:
      integrated = self.integrate(experiments, indexed)
    except Exception, e:
      print "Error integrating", tag, str(e)
      return

  def find_spots(self, datablock):
    from time import time
    from logging import info
    from dials.array_family import flex
    st = time()

    info('*' * 80)
    info('Finding Strong Spots')
    info('*' * 80)

    # Find the strong spots
    observed = flex.reflection_table.from_observations(datablock, self.params)

    # Save the reflections to file
    info('\n' + '-' * 80)
    if self.params.output.strong_filename:
      self.save_reflections(observed, self.params.output.strong_filename)

    info('')
    info('Time Taken = %f seconds' % (time() - st))
    return observed

  def index(self, datablock, reflections):
    from time import time
    from logging import info
    import copy
    st = time()

    info('*' * 80)
    info('Indexing Strong Spots')
    info('*' * 80)

    imagesets = datablock.extract_imagesets()

    params = copy.deepcopy(self.params)
    # don't do scan-varying refinement during indexing
    params.refinement.parameterisation.scan_varying = False

    from dials.algorithms.indexing.indexer import indexer_base
    idxr = indexer_base.from_parameters(
      reflections, imagesets,
      params=params)

    indexed = idxr.refined_reflections
    experiments = idxr.refined_experiments

    if self.params.output.indexed_filename:
      self.save_reflections(indexed, self.params.output.indexed_filename)

    info('')
    info('Time Taken = %f seconds' % (time() - st))
    return experiments, indexed

  def refine(self, experiments, centroids):
    print "Skipping refinement because the crystal orientation is refined during indexing"
    # Dump experiments to disk
    if self.params.output.refined_experiments_filename:
      from dxtbx.model.experiment.experiment_list import ExperimentListDumper
      dump = ExperimentListDumper(experiments)
      dump.as_json(self.params.output.refined_experiments_filename)

    return experiments

  def integrate(self, experiments, indexed):
    from time import time
    from logging import info

    st = time()

    info('*' * 80)
    info('Integrating Reflections')
    info('*' * 80)


    indexed,_ = self.process_reference(indexed)

    # Get the integrator from the input parameters
    info('Configuring integrator from input parameters')
    from dials.algorithms.profile_model.factory import ProfileModelFactory
    from dials.algorithms.integration.integrator import IntegratorFactory
    from dials.array_family import flex

    # Compute the profile model
    # Predict the reflections
    # Match the predictions with the reference
    # Create the integrator
    experiments = ProfileModelFactory.create(self.params, experiments, indexed)
    info("")
    info("=" * 80)
    info("")
    info("Predicting reflections")
    info("")
    predicted = flex.reflection_table.from_predictions_multi(
      experiments,
      dmin=self.params.prediction.d_min,
      dmax=self.params.prediction.d_max,
      margin=self.params.prediction.margin,
      force_static=self.params.prediction.force_static)
    predicted.match_with_reference(indexed)
    info("")
    integrator = IntegratorFactory.create(self.params, experiments, predicted)

    # Integrate the reflections
    integrated = integrator.integrate()

    # verify sigmas are sensible
    if 'intensity.prf.value' in integrated:
      if (integrated['intensity.prf.variance'] <= 0).count(True) > 0:
        raise Sorry("Found negative variances")
    if 'intensity.sum.value' in integrated:
      if (integrated['intensity.sum.variance'] <= 0).count(True) > 0:
        raise Sorry("Found negative variances")
      # apply detector gain to summation variances
      integrated['intensity.sum.variance'] *= self.params.integration.summation.detector_gain
    if 'background.sum.value' in integrated:
      if (integrated['background.sum.variance'] <= 0).count(True) > 0:
        raise Sorry("Found negative variances")
      # apply detector gain to background summation variances
      integrated['background.sum.variance'] *= self.params.integration.summation.detector_gain

    len_all = len(integrated)
    integrated = integrated.select(~integrated.get_flags(integrated.flags.foreground_includes_bad_pixels))
    print "Filtering %d reflections with at least one bad foreground pixel out of %d"%(len_all-len(integrated), len_all)

    if self.params.output.integrated_filename:
      # Save the reflections
      self.save_reflections(integrated, self.params.output.integrated_filename)

    self.write_integration_pickles(integrated, experiments)
    from dials.algorithms.indexing.stills_indexer import calc_2D_rmsd_and_displacements
    rmsd_indexed, _ = calc_2D_rmsd_and_displacements(indexed)
    rmsd_integrated, _ = calc_2D_rmsd_and_displacements(integrated)
    crystal_model = experiments.crystals()[0]
    print "Integrated. RMSD indexed,", rmsd_indexed, "RMSD integrated", rmsd_integrated, \
      "Final ML model: domain size angstroms: %f, half mosaicity degrees: %f"%(crystal_model._ML_domain_size_ang, crystal_model._ML_half_mosaicity_deg)

    info('')
    info('Time Taken = %f seconds' % (time() - st))
    return integrated

  def write_integration_pickles(self, integrated, experiments, callback = None):
    """
    Write a serialized python dictionary with integrated intensities and other information
    suitible for use by cxi.merge or prime.postrefine.
    @param integrated Reflection table with integrated intensities
    @param experiments Experiment list. One integration pickle for each experiment will be created.
    @param callback Deriving classes can use callback to make further modifications to the dictionary
    before it is serialized. Callback should be a function with this signature:
    def functionname(params, outfile, frame), where params is the phil scope, outfile is the path
    to the pickle that will be saved, and frame is the python dictionary to be serialized.
    """
    try:
      picklefilename = self.params.output.integration_pickle
    except AttributeError:
      return

    if self.params.output.integration_pickle is not None:

      from libtbx import easy_pickle
      import os
      from xfel.command_line.frame_extractor import ConstructFrame
      from dials.array_family import flex

      # Split everything into separate experiments for pickling
      for e_number in xrange(len(experiments)):
        experiment = experiments[e_number]
        e_selection = integrated['id'] == e_number
        reflections = integrated.select(e_selection)

        frame = ConstructFrame(reflections, experiment).make_frame()
        frame["pixel_size"] = experiment.detector[0].get_pixel_size()[0]

        try:
          # if the data was a file on disc, get the path
          event_timestamp = os.path.splitext(experiments[0].imageset.paths()[0])[0]
        except NotImplementedError:
          # if the data is in memory only, check if the reader set a timestamp on the format object
          event_timestamp = experiment.imageset.reader().get_format(0).timestamp
        event_timestamp = os.path.basename(event_timestamp)
        if event_timestamp.find("shot-")==0:
           event_timestamp = os.path.splitext(event_timestamp)[0] # micromanage the file name
        if hasattr(self.params.output, "output_dir"):
          outfile = os.path.join(self.params.output.output_dir, self.params.output.integration_pickle%(e_number,event_timestamp))
        else:
          outfile = os.path.join(os.path.dirname(self.params.output.integration_pickle), self.params.output.integration_pickle%(e_number,event_timestamp))

        if callback is not None:
          callback(self.params, outfile, frame)

        easy_pickle.dump(outfile, frame)

  def process_reference(self, reference):
    ''' Load the reference spots. '''
    from dials.array_family import flex
    from logging import info
    from time import time
    if reference is None:
      return None, None
    st = time()
    assert("miller_index" in reference)
    assert("id" in reference)
    info('Processing reference reflections')
    info(' read %d strong spots' % len(reference))
    mask = reference.get_flags(reference.flags.indexed)
    rubbish = reference.select(mask == False)
    if mask.count(False) > 0:
      reference.del_selected(mask == False)
      info(' removing %d unindexed reflections' %  mask.count(True))
    if len(reference) == 0:
      raise Sorry('''
        Invalid input for reference reflections.
        Expected > %d indexed spots, got %d
      ''' % (0, len(reference)))
    mask = reference['miller_index'] == (0, 0, 0)
    if mask.count(True) > 0:
      rubbish.extend(reference.select(mask))
      reference.del_selected(mask)
      info(' removing %d reflections with hkl (0,0,0)' %  mask.count(True))
    mask = reference['id'] < 0
    if mask.count(True) > 0:
      raise Sorry('''
        Invalid input for reference reflections.
        %d reference spots have an invalid experiment id
      ''' % mask.count(True))
    info(' using %d indexed reflections' % len(reference))
    info(' found %d junk reflections' % len(rubbish))
    info(' time taken: %g' % (time() - st))
    return reference, rubbish

  def save_reflections(self, reflections, filename):
    ''' Save the reflections to file. '''
    from logging import info
    from time import time
    st = time()
    info('Saving %d reflections to %s' % (len(reflections), filename))
    reflections.as_pickle(filename)
    info(' time taken: %g' % (time() - st))

if __name__ == '__main__':
  from dials.util import halraiser
  try:
    script = Script()
    script.run()
  except Exception as e:
    halraiser(e)

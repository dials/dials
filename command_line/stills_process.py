#!/usr/bin/env python
#
# LIBTBX_SET_DISPATCHER_NAME dials.stills_process

from __future__ import absolute_import, division
import libtbx.load_env
import logging
logger = logging.getLogger(libtbx.env.dispatcher_name)

from libtbx.utils import Abort, Sorry
from dxtbx.datablock import DataBlockFactory
import os

help_message = '''
DIALS script for processing still images. Import, index, refine, and integrate are all done for each image
seperately.
'''

from libtbx.phil import parse
control_phil_str = '''
  verbosity = 1
    .type = int(value_min=0)
    .help = "The verbosity level"

  dispatch {
    pre_import = False
      .type = bool
      .expert_level = 2
      .help = If True, before processing import all the data. Needed only if processing \
              multiple multi-image files at once (not a recommended use case)
    refine = False
      .expert_level = 2
      .type = bool
      .help = If True, after indexing, refine the experimental models
    squash_errors = True
      .expert_level = 2
      .type = bool
      .help = If True, if an image fails to process, continue to the next image. \
              otherwise, halt processing and show the error.
  }

  output {
    output_dir = .
      .type = str
      .help = Directory output files will be placed
    composite_output = False
      .type = bool
      .help = If True, save one set of json/pickle files per process, where each is a \
              concatenated list of all the successful events examined by that process. \
              If False, output a separate json/pickle file per image (generates a \
              lot of files).
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
    integrated_experiments_filename = %s_integrated_experiments.json
      .type = str
      .help = The filename for saving final experimental models.
    profile_filename = None
      .type = str
      .help = The filename for output reflection profile parameters
    integration_pickle = int-%d-%s.pickle
      .type = str
      .help = Filename for cctbx.xfel-style integration pickle files
  }

  mp {
    method = *multiprocessing sge lsf pbs mpi
      .type = choice
      .help = "The multiprocessing method to use"
    nproc = 1
      .type = int(value_min=1)
      .help = "The number of processes to use."
  }
'''

dials_phil_str = '''
  input {
    reference_geometry = None
      .type = str
      .help = Provide an experiments.json file with exactly one detector model. Data processing will use \
              that geometry instead of the geometry found in the image headers.
  }

  output {
    shoeboxes = True
      .type = bool
      .help = Save the raw pixel values inside the reflection shoeboxes during spotfinding.
  }

  include scope dials.util.options.geometry_phil_scope
  include scope dials.algorithms.spot_finding.factory.phil_scope
  include scope dials.algorithms.indexing.indexer.index_only_phil_scope
  include scope dials.algorithms.refinement.refiner.phil_scope
  include scope dials.algorithms.integration.integrator.phil_scope
  include scope dials.algorithms.profile_model.factory.phil_scope
  include scope dials.algorithms.spot_prediction.reflection_predictor.phil_scope
  include scope dials.algorithms.integration.stills_significance_filter.phil_scope

  indexing {
    stills {
      method_list = None
        .type = strings
        .help = List of indexing methods. If indexing fails with first method, indexing will be \
                attempted with the next, and so forth
    }
  }

  integration {
    include scope dials.algorithms.integration.kapton_correction.absorption_phil_scope
  }
'''

program_defaults_phil_str = '''
indexing {
  method = fft1d
}
refinement {
  parameterisation {
    auto_reduction {
      min_nref_per_parameter = 1
      action = fix
    }
    beam.fix = all
    detector.fix = all
  }
  reflections {
    weighting_strategy.override = stills
    outlier.algorithm = null
  }
}
integration {
  integrator = stills
  profile.fitting = False
  background {
    algorithm = simple
    simple {
      outlier.algorithm = plane
      model.algorithm = linear2d
    }
  }
}
profile.gaussian_rs.min_spots.overall = 0
'''

phil_scope = parse(control_phil_str + dials_phil_str, process_includes=True).fetch(parse(program_defaults_phil_str))

def do_import(filename):
  logger.info("Loading %s"%os.path.basename(filename))
  datablocks = DataBlockFactory.from_filenames([filename])
  if len(datablocks) == 0:
    try:
      datablocks = DataBlockFactory.from_json_file(filename)
    except ValueError:
      raise Abort("Could not load %s"%filename)

  if len(datablocks) == 0:
    raise Abort("Could not load %s"%filename)
  if len(datablocks) > 1:
    raise Abort("Got multiple datablocks from file %s"%filename)

  # Ensure the indexer and downstream applications treat this as set of stills
  from dxtbx.imageset import ImageSet
  reset_sets = []

  for imageset in datablocks[0].extract_imagesets():
    imageset = ImageSet(imageset.data(), imageset.indices())
    imageset.set_scan(None)
    imageset.set_goniometer(None)
    reset_sets.append(imageset)

  return DataBlockFactory.from_imageset(reset_sets)[0]

class Script(object):
  '''A class for running the script.'''

  def __init__(self):
    '''Initialise the script.'''
    from dials.util.options import OptionParser
    import libtbx.load_env

    # The script usage
    usage = "usage: %s [options] [param.phil] filenames" % libtbx.env.dispatcher_name

    self.tag = None
    self.reference_detector = None

    # Create the parser
    self.parser = OptionParser(
      usage=usage,
      phil=phil_scope,
      epilog=help_message
      )

  def load_reference_geometry(self):
    if self.params.input.reference_geometry is None: return

    try:
      ref_datablocks = DataBlockFactory.from_json_file(self.params.input.reference_geometry, check_format=False)
    except Exception:
      ref_datablocks = None
    if ref_datablocks is None:
      from dxtbx.model.experiment_list import ExperimentListFactory
      try:
        ref_experiments = ExperimentListFactory.from_json_file(self.params.input.reference_geometry, check_format=False)
      except Exception:
        try:
          import dxtbx
          img = dxtbx.load(self.params.input.reference_geometry)
        except Exception:
          raise Sorry("Couldn't load geometry file %s"%self.params.input.reference_geometry)
        else:
          self.reference_detector = img.get_detector()
      else:
        assert len(ref_experiments.detectors()) == 1
        self.reference_detector = ref_experiments.detectors()[0]
    else:
      assert len(ref_datablocks) == 1 and len(ref_datablocks[0].unique_detectors()) == 1
      self.reference_detector = ref_datablocks[0].unique_detectors()[0]

  def run(self):
    '''Execute the script.'''
    from dials.util import log
    from time import time
    from libtbx import easy_mp
    import copy

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
      logger.info('The following parameters have been modified:\n')
      logger.info(diff_phil)

    for abs_params in self.params.integration.absorption_correction:
      if abs_params.apply:
        if not (self.params.integration.debug.output and not self.params.integration.debug.separate_files):
          raise Sorry('Shoeboxes must be saved to integration intermediates to apply an absorption correction. '\
            +'Set integration.debug.output=True, integration.debug.separate_files=False and '\
            +'integration.debug.delete_shoeboxes=True to temporarily store shoeboxes.')

    self.load_reference_geometry()
    from dials.command_line.dials_import import ManualGeometryUpdater
    update_geometry = ManualGeometryUpdater(params)

    # Import stuff
    logger.info("Loading files...")
    pre_import = params.dispatch.pre_import or len(all_paths) == 1
    if pre_import:
      # Handle still imagesets by breaking them apart into multiple datablocks
      # Further handle single file still imagesets (like HDF5) by tagging each
      # frame using its index

      datablocks = [do_import(path) for path in all_paths]
      if self.reference_detector is not None:
        from dxtbx.model import Detector
        for datablock in datablocks:
          for imageset in datablock.extract_imagesets():
            for i in range(len(imageset)):
              imageset.set_detector(
                Detector.from_dict(self.reference_detector.to_dict()),
                index=i)

      for datablock in datablocks:
        for imageset in datablock.extract_imagesets():
          update_geometry(imageset)

      indices = []
      basenames = []
      split_datablocks = []
      for datablock in datablocks:
        for imageset in datablock.extract_imagesets():
          paths = imageset.paths()
          for i in xrange(len(imageset)):
            subset = imageset[i:i+1]
            split_datablocks.append(DataBlockFactory.from_imageset(subset)[0])
            indices.append(i)
            basenames.append(os.path.splitext(os.path.basename(paths[i]))[0])
      tags = []
      for i, basename in zip(indices, basenames):
        if basenames.count(basename) > 1:
          tags.append("%s_%05d"%(basename, i))
        else:
          tags.append(basename)

      # Wrapper function
      def do_work(i, item_list):
        processor = Processor(copy.deepcopy(params), composite_tag = "%04d"%i)
        for item in item_list:
          processor.process_datablock(item[0], item[1])
        processor.finalize()

      iterable = zip(tags, split_datablocks)

    else:
      basenames = [os.path.splitext(os.path.basename(filename))[0] for filename in all_paths]
      tags = []
      for i, basename in enumerate(basenames):
        if basenames.count(basename) > 1:
          tags.append("%s_%05d"%(basename, i))
        else:
          tags.append(basename)

      # Wrapper function
      def do_work(i, item_list):
        processor = Processor(copy.deepcopy(params), composite_tag = "%04d"%i)
        for item in item_list:
          tag, filename = item

          datablock = do_import(filename)
          imagesets = datablock.extract_imagesets()
          if len(imagesets) == 0 or len(imagesets[0]) == 0:
            logger.info("Zero length imageset in file: %s"%filename)
            return
          if len(imagesets) > 1:
            raise Abort("Found more than one imageset in file: %s"%filename)
          if len(imagesets[0]) > 1:
            raise Abort("Found a multi-image file. Run again with pre_import=True")

          if self.reference_detector is not None:
            from dxtbx.model import Detector
            imagesets[0].set_detector(Detector.from_dict(self.reference_detector.to_dict()))

          update_geometry(imagesets[0])

          processor.process_datablock(tag, datablock)
        processor.finalize()

      iterable = zip(tags, all_paths)

    # Process the data
    if params.mp.method == 'mpi':
      from mpi4py import MPI
      comm = MPI.COMM_WORLD
      rank = comm.Get_rank() # each process in MPI has a unique id, 0-indexed
      size = comm.Get_size() # size: number of processes running in this job

      subset = [item for i, item in enumerate(iterable) if (i+rank)%size == 0]
      do_work((rank, subset))
    else:
      from dxtbx.command_line.image_average import splitit
      result = list(easy_mp.multi_core_run(
        myfunction=do_work,
        argstuples=list(enumerate(splitit(iterable, params.mp.nproc))),
        nproc=params.mp.nproc))
      error_list = [r[2] for r in result]
      if error_list.count(None) != len(error_list):
        print "Some processes failed excecution. Not all images may have processed. Error messages:"
        for error in error_list:
          if error is None: continue
          print error

    # Total Time
    logger.info("")
    logger.info("Total Time Taken = %f seconds" % (time() - st))

class Processor(object):
  def __init__(self, params, composite_tag = None):
    self.params = params
    self.composite_tag = composite_tag

    # The convention is to put %s in the phil parameter to add a tag to
    # each output datafile. Save the initial templates here.
    self.datablock_filename_template              = params.output.datablock_filename
    self.strong_filename_template                 = params.output.strong_filename
    self.indexed_filename_template                = params.output.indexed_filename
    self.refined_experiments_filename_template    = params.output.refined_experiments_filename
    self.integrated_filename_template             = params.output.integrated_filename
    self.integrated_experiments_filename_template = params.output.integrated_experiments_filename

    if params.output.composite_output:
      assert composite_tag is not None
      from dxtbx.model.experiment_list import ExperimentList
      from dials.array_family import flex
      #self.all_strong_reflections = flex.reflection_table() # no composite strong pickles yet
      self.all_indexed_experiments = ExperimentList()
      self.all_indexed_reflections = flex.reflection_table()
      self.all_integrated_experiments = ExperimentList()
      self.all_integrated_reflections = flex.reflection_table()
      self.all_int_pickle_filenames = []
      self.all_int_pickles = []

      self.setup_filenames(composite_tag)

  def setup_filenames(self, tag):
    # before processing, set output paths according to the templates
    if self.datablock_filename_template is not None and "%s" in self.datablock_filename_template:
      self.params.output.datablock_filename = os.path.join(self.params.output.output_dir, self.datablock_filename_template%("idx-" + tag))
    if self.strong_filename_template is not None and "%s" in self.strong_filename_template:
      self.params.output.strong_filename = os.path.join(self.params.output.output_dir, self.strong_filename_template%("idx-" + tag))
    if self.indexed_filename_template is not None and "%s" in self.indexed_filename_template:
      self.params.output.indexed_filename = os.path.join(self.params.output.output_dir, self.indexed_filename_template%("idx-" + tag))
    if self.refined_experiments_filename_template is not None and "%s" in self.refined_experiments_filename_template:
      self.params.output.refined_experiments_filename = os.path.join(self.params.output.output_dir, self.refined_experiments_filename_template%("idx-" + tag))
    if self.integrated_filename_template is not None and "%s" in self.integrated_filename_template:
      self.params.output.integrated_filename = os.path.join(self.params.output.output_dir, self.integrated_filename_template%("idx-" + tag))
    if self.integrated_experiments_filename_template is not None and "%s" in self.integrated_experiments_filename_template:
      self.params.output.integrated_experiments_filename = os.path.join(self.params.output.output_dir, self.integrated_experiments_filename_template%("idx-" + tag))

  def process_datablock(self, tag, datablock):
    import os

    if not self.params.output.composite_output:
      self.setup_filenames(tag)
    self.tag = tag

    if self.params.output.datablock_filename:
      from dxtbx.datablock import DataBlockDumper
      dump = DataBlockDumper(datablock)
      dump.as_json(self.params.output.datablock_filename)

    # Do the processing
    try:
      observed = self.find_spots(datablock)
    except Exception as e:
      print "Error spotfinding", tag, str(e)
      if not self.params.dispatch.squash_errors: raise
      return
    try:
      experiments, indexed = self.index(datablock, observed)
    except Exception as e:
      print "Couldn't index", tag, str(e)
      if not self.params.dispatch.squash_errors: raise
      return
    try:
      experiments, indexed = self.refine(experiments, indexed)
    except Exception as e:
      print "Error refining", tag, str(e)
      if not self.params.dispatch.squash_errors: raise
      return
    try:
      integrated = self.integrate(experiments, indexed)
    except Exception as e:
      print "Error integrating", tag, str(e)
      if not self.params.dispatch.squash_errors: raise
      return

  def find_spots(self, datablock):
    from time import time
    from dials.array_family import flex
    st = time()

    logger.info('*' * 80)
    logger.info('Finding Strong Spots')
    logger.info('*' * 80)

    # Find the strong spots
    observed = flex.reflection_table.from_observations(datablock, self.params)

    # Reset z coordinates for dials.image_viewer; see Issues #226 for details
    xyzobs = observed['xyzobs.px.value']
    for i in xrange(len(xyzobs)):
      xyzobs[i] = (xyzobs[i][0], xyzobs[i][1], 0)
    bbox = observed['bbox']
    for i in xrange(len(bbox)):
      bbox[i] = (bbox[i][0], bbox[i][1], bbox[i][2], bbox[i][3], 0, 1)

    if self.params.output.composite_output:
      pass # no composite strong pickles yet
    else:
      # Save the reflections to file
      logger.info('\n' + '-' * 80)
      if self.params.output.strong_filename:
        self.save_reflections(observed, self.params.output.strong_filename)

    logger.info('')
    logger.info('Time Taken = %f seconds' % (time() - st))
    return observed

  def index(self, datablock, reflections):
    from dials.algorithms.indexing.indexer import indexer_base
    from time import time
    import copy
    st = time()

    logger.info('*' * 80)
    logger.info('Indexing Strong Spots')
    logger.info('*' * 80)

    imagesets = datablock.extract_imagesets()

    params = copy.deepcopy(self.params)
    # don't do scan-varying refinement during indexing
    params.refinement.parameterisation.scan_varying = False

    if hasattr(self, 'known_crystal_models'):
      known_crystal_models = self.known_crystal_models
    else:
      known_crystal_models = None

    if params.indexing.stills.method_list is None:
      idxr = indexer_base.from_parameters(
        reflections, imagesets, known_crystal_models=known_crystal_models,
        params=params)
      idxr.index()
    else:
      indexing_error = None
      for method in params.indexing.stills.method_list:
        params.indexing.method = method
        try:
          idxr = indexer_base.from_parameters(
            reflections, imagesets,
            params=params)
          idxr.index()
        except Exception as e:
          logger.info("Couldn't index using method %s"%method)
          if indexing_error is None:
            if e is None:
              e = Exception("Couldn't index using method %s"%method)
            indexing_error = e
        else:
          indexing_error = None
          break
      if indexing_error is not None:
        raise indexing_error

    indexed = idxr.refined_reflections
    experiments = idxr.refined_experiments

    if known_crystal_models is not None:
      from dials.array_family import flex
      filtered = flex.reflection_table()
      for idx in set(indexed['miller_index']):
        sel = indexed['miller_index'] == idx
        if sel.count(True) == 1:
          filtered.extend(indexed.select(sel))
      logger.info("Filtered duplicate reflections, %d out of %d remaining"%(len(filtered),len(indexed)))
      print "Filtered duplicate reflections, %d out of %d remaining"%(len(filtered),len(indexed))
      indexed = filtered

    logger.info('')
    logger.info('Time Taken = %f seconds' % (time() - st))
    return experiments, indexed

  def refine(self, experiments, centroids):
    if self.params.dispatch.refine:
      from dials.algorithms.refinement import RefinerFactory
      from time import time
      st = time()

      logger.info('*' * 80)
      logger.info('Refining Model')
      logger.info('*' * 80)

      refiner = RefinerFactory.from_parameters_data_experiments(
        self.params, centroids, experiments)

      refiner.run()
      experiments = refiner.get_experiments()
      predicted = refiner.predict_for_indexed()
      centroids['xyzcal.mm'] = predicted['xyzcal.mm']
      centroids['entering'] = predicted['entering']
      centroids = centroids.select(refiner.selection_used_for_refinement())

      # Re-estimate mosaic estimates
      from dials.algorithms.indexing.nave_parameters import nave_parameters
      nv = nave_parameters(params = self.params, experiments=experiments, reflections=centroids, refinery=refiner, graph_verbose=False)
      nv()
      acceptance_flags_nv = nv.nv_acceptance_flags
      centroids = centroids.select(acceptance_flags_nv)

    if self.params.output.composite_output:
      if self.params.output.refined_experiments_filename or self.params.output.indexed_filename:
        assert self.params.output.refined_experiments_filename is not None and self.params.output.indexed_filename is not None
        from dials.array_family import flex
        n = len(self.all_indexed_experiments)
        self.all_indexed_experiments.extend(experiments)
        for i, experiment in enumerate(experiments):
          refls = centroids.select(centroids['id'] == i)
          refls['id'] = flex.int(len(refls), n)
          self.all_indexed_reflections.extend(refls)
          n += 1
    else:
      # Dump experiments to disk
      if self.params.output.refined_experiments_filename:
        from dxtbx.model.experiment_list import ExperimentListDumper
        dump = ExperimentListDumper(experiments)
        dump.as_json(self.params.output.refined_experiments_filename)

      if self.params.output.indexed_filename:
        self.save_reflections(centroids, self.params.output.indexed_filename)

    if self.params.dispatch.refine:
      logger.info('')
      logger.info('Time Taken = %f seconds' % (time() - st))

    return experiments, centroids

  def integrate(self, experiments, indexed):
    from time import time

    st = time()

    logger.info('*' * 80)
    logger.info('Integrating Reflections')
    logger.info('*' * 80)


    indexed,_ = self.process_reference(indexed)

    # Get the integrator from the input parameters
    logger.info('Configuring integrator from input parameters')
    from dials.algorithms.profile_model.factory import ProfileModelFactory
    from dials.algorithms.integration.integrator import IntegratorFactory
    from dials.array_family import flex

    # Compute the profile model
    # Predict the reflections
    # Match the predictions with the reference
    # Create the integrator
    experiments = ProfileModelFactory.create(self.params, experiments, indexed)
    logger.info("")
    logger.info("=" * 80)
    logger.info("")
    logger.info("Predicting reflections")
    logger.info("")
    predicted = flex.reflection_table.from_predictions_multi(
      experiments,
      dmin=self.params.prediction.d_min,
      dmax=self.params.prediction.d_max,
      margin=self.params.prediction.margin,
      force_static=self.params.prediction.force_static)
    predicted.match_with_reference(indexed)
    logger.info("")
    integrator = IntegratorFactory.create(self.params, experiments, predicted)

    # Integrate the reflections
    integrated = integrator.integrate()

    # correct integrated intensities for absorption correction, if necessary
    for abs_params in self.params.integration.absorption_correction:
      if abs_params.apply and abs_params.algorithm == "fuller_kapton":
        from dials.algorithms.integration.kapton_correction import multi_kapton_correction
        experiments, integrated = multi_kapton_correction(experiments, integrated,
          abs_params.fuller_kapton, logger=logger)()

    if self.params.significance_filter.enable:
      from dials.algorithms.integration.stills_significance_filter import SignificanceFilter
      sig_filter = SignificanceFilter(self.params)
      refls = sig_filter(experiments, integrated)
      logger.info("Removed %d reflections out of %d when applying significance filter"%(len(integrated)-len(refls), len(integrated)))
      if len(refls) == 0:
        raise Sorry("No reflections left after applying significance filter")
      integrated = refls

    # Delete the shoeboxes used for intermediate calculations, if requested
    if self.params.integration.debug.delete_shoeboxes and 'shoebox' in integrated:
      del integrated['shoebox']

    if self.params.output.composite_output:
      if self.params.output.integrated_experiments_filename or self.params.output.integrated_filename:
        assert self.params.output.integrated_experiments_filename is not None and self.params.output.integrated_filename is not None
        from dials.array_family import flex
        n = len(self.all_integrated_experiments)
        self.all_integrated_experiments.extend(experiments)
        for i, experiment in enumerate(experiments):
          refls = integrated.select(integrated['id'] == i)
          refls['id'] = flex.int(len(refls), n)
          self.all_integrated_reflections.extend(refls)
          n += 1
    else:
      # Dump experiments to disk
      if self.params.output.integrated_experiments_filename:
        from dxtbx.model.experiment_list import ExperimentListDumper
        dump = ExperimentListDumper(experiments)
        dump.as_json(self.params.output.integrated_experiments_filename)

      if self.params.output.integrated_filename:
        # Save the reflections
        self.save_reflections(integrated, self.params.output.integrated_filename)

    self.write_integration_pickles(integrated, experiments)
    from dials.algorithms.indexing.stills_indexer import calc_2D_rmsd_and_displacements

    rmsd_indexed, _ = calc_2D_rmsd_and_displacements(indexed)
    log_str = "RMSD indexed (px): %f\n"%(rmsd_indexed)
    for i in xrange(6):
      bright_integrated = integrated.select((integrated['intensity.sum.value']/flex.sqrt(integrated['intensity.sum.variance']))>=i)
      if len(bright_integrated) > 0:
        rmsd_integrated, _ = calc_2D_rmsd_and_displacements(bright_integrated)
      else:
        rmsd_integrated = 0
      log_str += "N reflections integrated at I/sigI >= %d: % 4d, RMSD (px): %f\n"%(i, len(bright_integrated), rmsd_integrated)

    for crystal_model in experiments.crystals():
      if hasattr(crystal_model, 'get_domain_size_ang'):
        log_str += ". Final ML model: domain size angstroms: %f, half mosaicity degrees: %f"%(crystal_model.get_domain_size_ang(), crystal_model.get_half_mosaicity_deg())

    logger.info(log_str)

    logger.info('')
    logger.info('Time Taken = %f seconds' % (time() - st))
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

        if not hasattr(self, 'tag') or self.tag is None:
          try:
            # if the data was a file on disc, get the path
            event_timestamp = os.path.splitext(experiments[0].imageset.paths()[0])[0]
          except NotImplementedError:
            # if the data is in memory only, check if the reader set a timestamp on the format object
            event_timestamp = experiment.imageset.reader().get_format(0).timestamp
          event_timestamp = os.path.basename(event_timestamp)
          if event_timestamp.find("shot-")==0:
             event_timestamp = os.path.splitext(event_timestamp)[0] # micromanage the file name
        else:
          event_timestamp = self.tag
        if hasattr(self.params.output, "output_dir"):
          outfile = os.path.join(self.params.output.output_dir, self.params.output.integration_pickle%(e_number,event_timestamp))
        else:
          outfile = os.path.join(os.path.dirname(self.params.output.integration_pickle), self.params.output.integration_pickle%(e_number,event_timestamp))

        if callback is not None:
          callback(self.params, outfile, frame)

        if self.params.output.composite_output:
          self.all_int_pickle_filenames.append(os.path.basename(outfile))
          self.all_int_pickles.append(frame)
        else:
          easy_pickle.dump(outfile, frame)

  def process_reference(self, reference):
    ''' Load the reference spots. '''
    from dials.array_family import flex
    from time import time
    if reference is None:
      return None, None
    st = time()
    assert("miller_index" in reference)
    assert("id" in reference)
    logger.info('Processing reference reflections')
    logger.info(' read %d strong spots' % len(reference))
    mask = reference.get_flags(reference.flags.indexed)
    rubbish = reference.select(mask == False)
    if mask.count(False) > 0:
      reference.del_selected(mask == False)
      logger.info(' removing %d unindexed reflections' %  mask.count(True))
    if len(reference) == 0:
      raise Sorry('''
        Invalid input for reference reflections.
        Expected > %d indexed spots, got %d
      ''' % (0, len(reference)))
    mask = reference['miller_index'] == (0, 0, 0)
    if mask.count(True) > 0:
      rubbish.extend(reference.select(mask))
      reference.del_selected(mask)
      logger.info(' removing %d reflections with hkl (0,0,0)' %  mask.count(True))
    mask = reference['id'] < 0
    if mask.count(True) > 0:
      raise Sorry('''
        Invalid input for reference reflections.
        %d reference spots have an invalid experiment id
      ''' % mask.count(True))
    logger.info(' using %d indexed reflections' % len(reference))
    logger.info(' found %d junk reflections' % len(rubbish))
    logger.info(' time taken: %g' % (time() - st))
    return reference, rubbish

  def save_reflections(self, reflections, filename):
    ''' Save the reflections to file. '''
    from time import time
    st = time()
    logger.info('Saving %d reflections to %s' % (len(reflections), filename))
    reflections.as_pickle(filename)
    logger.info(' time taken: %g' % (time() - st))

  def finalize(self):
    ''' Perform any final operations '''
    if self.params.output.composite_output:
      # Dump composite files to disk
      if len(self.all_indexed_experiments) > 0 and self.params.output.refined_experiments_filename:
        from dxtbx.model.experiment_list import ExperimentListDumper
        dump = ExperimentListDumper(self.all_indexed_experiments)
        dump.as_json(self.params.output.refined_experiments_filename)

      if len(self.all_indexed_reflections) > 0 and self.params.output.indexed_filename:
        self.save_reflections(self.all_indexed_reflections, self.params.output.indexed_filename)

      if len(self.all_integrated_experiments) > 0 and self.params.output.integrated_experiments_filename:
        from dxtbx.model.experiment_list import ExperimentListDumper
        dump = ExperimentListDumper(self.all_integrated_experiments)
        dump.as_json(self.params.output.integrated_experiments_filename)

      if len(self.all_integrated_reflections) > 0 and self.params.output.integrated_filename:
        self.save_reflections(self.all_integrated_reflections, self.params.output.integrated_filename)

      # Create a tar archive of the integration dictionary pickles
      if len(self.all_int_pickles) > 0 and self.params.output.integration_pickle:
        import tarfile, StringIO, time, cPickle as pickle
        tar_template_integration_pickle = self.params.output.integration_pickle.replace('%d', '%s')
        outfile = os.path.join(self.params.output.output_dir, tar_template_integration_pickle%('x',self.composite_tag)) + ".tar"
        tar = tarfile.TarFile(outfile,"w")
        for i, (fname, d) in enumerate(zip(self.all_int_pickle_filenames, self.all_int_pickles)):
          string = StringIO.StringIO(pickle.dumps(d, protocol=2))
          info = tarfile.TarInfo(name=fname)
          info.size=len(string.buf)
          info.mtime = time.time()
          tar.addfile(tarinfo=info, fileobj=string)
        tar.close()

if __name__ == '__main__':
  from dials.util import halraiser
  try:
    script = Script()
    script.run()
  except Exception as e:
    halraiser(e)

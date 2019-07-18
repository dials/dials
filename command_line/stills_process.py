#!/usr/bin/env python

from __future__ import absolute_import, division, print_function

import copy
import logging
import os
import sys
import tarfile
import time
import six.moves.cPickle as pickle
from six.moves import StringIO

import dials.util
from dials.array_family import flex
from dxtbx.model.experiment_list import ExperimentListFactory
from dxtbx.model.experiment_list import ExperimentList
from libtbx.utils import Abort, Sorry

logger = logging.getLogger("dials.command_line.stills_process")


help_message = """
DIALS script for processing still images. Import, index, refine, and integrate are all done for each image
seperately.
"""

from libtbx.phil import parse

control_phil_str = """
  verbosity = 0
    .type = int(value_min=0)
    .help = "The verbosity level"

  input {
    file_list = None
      .type = path
      .help = Path to a list of images
  }

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
    find_spots = True
      .expert_level = 2
      .type = bool
      .help = Whether to do spotfinding. Needed for indexing/integration
    index = True
      .expert_level = 2
      .type = bool
      .help = Attempt to index images. find_spots also needs to be True for this to work
    integrate = True
      .expert_level = 2
      .type = bool
      .help = Integrate indexed images. Ignored if index=False or find_spots=False
    hit_finder{
      enable = True
        .type = bool
        .help = Whether to do hitfinding. hit_finder=False: process all images
      minimum_number_of_reflections = 16
        .type = int
        .help = If the number of strong reflections on an image is less than this, and \
                 the hitfinder is enabled, discard this image.
      maximum_number_of_reflections = None
       .type = int
       .help = If specified, ignores images with more than this many number of reflections
    }
  }

  output {
    output_dir = .
      .type = str
      .help = Directory output files will be placed
    composite_output = False
      .type = bool
      .help = If True, save one set of experiment/reflection files per process, where each is a \
              concatenated list of all the successful events examined by that process. \
              If False, output a separate experiment/reflection file per image (generates a \
              lot of files).
    logging_dir = None
      .type = str
      .help = Directory output log files will be placed
    experiments_filename = %s_imported.expt
      .type = str
      .help = The filename for output experiments
    strong_filename = %s_strong.refl
      .type = str
      .help = The filename for strong reflections from spot finder output.
    indexed_filename = %s_indexed.refl
      .type = str
      .help = The filename for indexed reflections.
    refined_experiments_filename = %s_refined.expt
      .type = str
      .help = The filename for saving refined experimental models
    integrated_filename = %s_integrated.refl
      .type = str
      .help = The filename for final integrated reflections.
    integrated_experiments_filename = %s_integrated.expt
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
    glob = None
      .type = str
      .help = For MPI, for multifile data, mandatory blobs giving file paths
      .multiple = True
    composite_stride = None
      .type = int
      .help = For MPI, if using composite mode, specify how many ranks to    \
              aggregate data from.  For example, if you have 100 processes,  \
              composite mode will output N*100 files, where N is the number  \
              of file types (expt, refl, etc). If you specify stride = 25, \
              then each group of 25 process will send their results to 4     \
              processes and only N*4 files will be created. Ideally, match   \
              stride to the number of processors per node.
  }
"""

dials_phil_str = """
  input {
    reference_geometry = None
      .type = str
      .help = Provide an models.expt file with exactly one detector model. Data processing will use \
              that geometry instead of the geometry found in the image headers.
  }

  output {
    shoeboxes = True
      .type = bool
      .help = Save the raw pixel values inside the reflection shoeboxes during spotfinding.
  }

  include scope dials.util.options.geometry_phil_scope
  include scope dials.algorithms.spot_finding.factory.phil_scope
  include scope dials.algorithms.indexing.indexer.phil_scope
  indexing {
      include scope dials.algorithms.indexing.lattice_search.basis_vector_search_phil_scope
  }
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
"""

program_defaults_phil_str = """
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
"""

phil_scope = parse(control_phil_str + dials_phil_str, process_includes=True).fetch(
    parse(program_defaults_phil_str)
)


def do_import(filename, load_models=True):
    logger.info("Loading %s" % os.path.basename(filename))
    experiments = ExperimentListFactory.from_filenames([filename], load_models=False)
    if len(experiments) == 0:
        try:
            experiments = ExperimentListFactory.from_json_file(filename)
        except ValueError:
            raise Abort("Could not load %s" % filename)

    if len(experiments) == 0:
        raise Abort("Could not load %s" % filename)

    from dxtbx.imageset import ImageSetFactory

    for experiment in experiments:
        if load_models:
            experiment.load_models()
        imageset = ImageSetFactory.imageset_from_anyset(experiment.imageset)
        imageset.set_scan(None)
        imageset.set_goniometer(None)
        experiment.imageset = imageset
        experiment.scan = None
        experiment.goniometer = None

    return experiments


class Script(object):
    """A class for running the script."""

    def __init__(self):
        """Initialise the script."""
        from dials.util.options import OptionParser

        # The script usage
        usage = "usage: dials.stills_process [options] [param.phil] filenames"

        self.tag = None
        self.reference_detector = None

        # Create the parser
        self.parser = OptionParser(usage=usage, phil=phil_scope, epilog=help_message)

    def load_reference_geometry(self):
        if self.params.input.reference_geometry is None:
            return

        from dxtbx.model.experiment_list import ExperimentListFactory

        try:
            ref_experiments = ExperimentListFactory.from_json_file(
                self.params.input.reference_geometry, check_format=False
            )
        except Exception:
            try:
                import dxtbx

                img = dxtbx.load(self.params.input.reference_geometry)
            except Exception:
                raise Sorry(
                    "Couldn't load geometry file %s"
                    % self.params.input.reference_geometry
                )
            else:
                self.reference_detector = img.get_detector()
        else:
            assert len(ref_experiments.detectors()) == 1
            self.reference_detector = ref_experiments.detectors()[0]

    def run(self):
        """Execute the script."""
        from dials.util import log
        from libtbx import easy_mp

        # Parse the command line
        params, options, all_paths = self.parser.parse_args(
            show_diff_phil=False, return_unhandled=True, quick_parse=True
        )

        if not all_paths and params.input.file_list is not None:
            all_paths.extend(
                [path.strip() for path in open(params.input.file_list).readlines()]
            )

        # Check we have some filenames
        if not all_paths:
            self.parser.print_help()
            return

        # Mask validation
        for mask_path in params.spotfinder.lookup.mask, params.integration.lookup.mask:
            if mask_path is not None and not os.path.isfile(mask_path):
                raise Sorry("Mask %s not found" % mask_path)

        # Save the options
        self.options = options
        self.params = params

        st = time.time()

        # Configure logging
        log.config(
            params.verbosity, info="dials.process.log", debug="dials.process.debug.log"
        )

        bad_phils = [f for f in all_paths if os.path.splitext(f)[1] == ".phil"]
        if len(bad_phils) > 0:
            self.parser.print_help()
            logger.error(
                "Error: the following phil files were not understood: %s"
                % (", ".join(bad_phils))
            )
            return

        # Log the diff phil
        diff_phil = self.parser.diff_phil.as_str()
        if diff_phil != "":
            logger.info("The following parameters have been modified:\n")
            logger.info(diff_phil)

        for abs_params in self.params.integration.absorption_correction:
            if abs_params.apply:
                if not (
                    self.params.integration.debug.output
                    and not self.params.integration.debug.separate_files
                ):
                    raise Sorry(
                        "Shoeboxes must be saved to integration intermediates to apply an absorption correction. "
                        + "Set integration.debug.output=True, integration.debug.separate_files=False and "
                        + "integration.debug.delete_shoeboxes=True to temporarily store shoeboxes."
                    )

        self.load_reference_geometry()
        from dials.command_line.dials_import import ManualGeometryUpdater

        update_geometry = ManualGeometryUpdater(params)

        # Import stuff
        logger.info("Loading files...")
        pre_import = params.dispatch.pre_import or len(all_paths) == 1
        if pre_import:
            # Handle still imagesets by breaking them apart into multiple experiments
            # Further handle single file still imagesets (like HDF5) by tagging each
            # frame using its index

            experiments = ExperimentList()
            for path in all_paths:
                experiments.extend(do_import(path, load_models=False))

            indices = []
            basenames = []
            split_experiments = []
            for i, imageset in enumerate(experiments.imagesets()):
                assert len(imageset) == 1
                paths = imageset.paths()
                indices.append(i)
                basenames.append(os.path.splitext(os.path.basename(paths[0]))[0])
                split_experiments.append(experiments[i : i + 1])
            tags = []
            for i, basename in zip(indices, basenames):
                if basenames.count(basename) > 1:
                    tags.append("%s_%05d" % (basename, i))
                else:
                    tags.append(basename)

            # Wrapper function
            def do_work(i, item_list):
                processor = Processor(
                    copy.deepcopy(params), composite_tag="%04d" % i, rank=i
                )

                for item in item_list:
                    try:
                        assert len(item[1]) == 1
                        experiment = item[1][0]
                        experiment.load_models()
                        imageset = experiment.imageset
                        update_geometry(imageset)
                        experiment.beam = imageset.get_beam()
                        experiment.detector = imageset.get_detector()
                    except RuntimeError as e:
                        logger.warning(
                            "Error updating geometry on item %s, %s"
                            % (str(item[0]), str(e))
                        )
                        continue

                    if self.reference_detector is not None:
                        from dxtbx.model import Detector

                        experiment = item[1][0]
                        imageset = experiment.imageset
                        imageset.set_detector(
                            Detector.from_dict(self.reference_detector.to_dict())
                        )
                        experiment.detector = imageset.get_detector()

                    processor.process_experiments(item[0], item[1])
                processor.finalize()

            iterable = list(zip(tags, split_experiments))

        else:
            basenames = [
                os.path.splitext(os.path.basename(filename))[0]
                for filename in all_paths
            ]
            tags = []
            for i, basename in enumerate(basenames):
                if basenames.count(basename) > 1:
                    tags.append("%s_%05d" % (basename, i))
                else:
                    tags.append(basename)

            # Wrapper function
            def do_work(i, item_list):
                processor = Processor(
                    copy.deepcopy(params), composite_tag="%04d" % i, rank=i
                )
                for item in item_list:
                    tag, filename = item

                    experiments = do_import(filename, load_models=True)
                    imagesets = experiments.imagesets()
                    if len(imagesets) == 0 or len(imagesets[0]) == 0:
                        logger.info("Zero length imageset in file: %s" % filename)
                        return
                    if len(imagesets) > 1:
                        raise Abort(
                            "Found more than one imageset in file: %s" % filename
                        )
                    if len(imagesets[0]) > 1:
                        raise Abort(
                            "Found a multi-image file. Run again with pre_import=True"
                        )

                    try:
                        update_geometry(imagesets[0])
                        experiment = experiments[0]
                        experiment.beam = imagesets[0].get_beam()
                        experiment.detector = imagesets[0].get_detector()
                    except RuntimeError as e:
                        logger.warning(
                            "Error updating geometry on item %s, %s" % (tag, str(e))
                        )
                        continue

                    if self.reference_detector is not None:
                        from dxtbx.model import Detector

                        imageset = experiments[0].imageset
                        imageset.set_detector(
                            Detector.from_dict(self.reference_detector.to_dict())
                        )
                        experiments[0].detector = imageset.get_detector()

                    processor.process_experiments(tag, experiments)
                processor.finalize()

            iterable = list(zip(tags, all_paths))

        # Process the data
        if params.mp.method == "mpi":
            from mpi4py import MPI

            comm = MPI.COMM_WORLD
            rank = comm.Get_rank()  # each process in MPI has a unique id, 0-indexed
            size = comm.Get_size()  # size: number of processes running in this job

            # Configure the logging
            if params.output.logging_dir is None:
                info_path = ""
                debug_path = ""
            else:
                log_path = os.path.join(
                    params.output.logging_dir, "log_rank%04d.out" % rank
                )
                error_path = os.path.join(
                    params.output.logging_dir, "error_rank%04d.out" % rank
                )
                print("Redirecting stdout to %s" % log_path)
                print("Redirecting stderr to %s" % error_path)
                sys.stdout = open(log_path, "a", buffering=0)
                sys.stderr = open(error_path, "a", buffering=0)
                print("Should be redirected now")

                info_path = os.path.join(
                    params.output.logging_dir, "info_rank%04d.out" % rank
                )
                debug_path = os.path.join(
                    params.output.logging_dir, "debug_rank%04d.out" % rank
                )

            from dials.util import log

            log.config(params.verbosity, info=info_path, debug=debug_path)

            if size <= 2:  # client/server only makes sense for n>2
                subset = [
                    item for i, item in enumerate(iterable) if (i + rank) % size == 0
                ]
                do_work(rank, subset)
            else:
                if rank == 0:
                    # server process
                    for item in iterable:
                        print("Getting next available process")
                        rankreq = comm.recv(source=MPI.ANY_SOURCE)
                        print("Process %s is ready, sending %s\n" % (rankreq, item[0]))
                        comm.send(item, dest=rankreq)
                    # send a stop command to each process
                    print("MPI DONE, sending stops\n")
                    for rankreq in range(size - 1):
                        rankreq = comm.recv(source=MPI.ANY_SOURCE)
                        print("Sending stop to %d\n" % rankreq)
                        comm.send("endrun", dest=rankreq)
                    print("All stops sent.")
                else:
                    # client process
                    while True:
                        # inform the server this process is ready for an event
                        print("Rank %d getting next task" % rank)
                        comm.send(rank, dest=0)
                        print("Rank %d waiting for response" % rank)
                        item = comm.recv(source=0)
                        if item == "endrun":
                            print("Rank %d received endrun" % rank)
                            break
                        print("Rank %d beginning processing" % rank)
                        try:
                            do_work(rank, [item])
                        except Exception as e:
                            print(
                                "Rank %d unhandled exception processing event" % rank,
                                str(e),
                            )
                        print("Rank %d event processed" % rank)
        else:
            from dxtbx.command_line.image_average import splitit

            if params.mp.nproc == 1:
                do_work(0, iterable)
            else:
                result = list(
                    easy_mp.multi_core_run(
                        myfunction=do_work,
                        argstuples=list(enumerate(splitit(iterable, params.mp.nproc))),
                        nproc=params.mp.nproc,
                    )
                )
                error_list = [r[2] for r in result]
                if error_list.count(None) != len(error_list):
                    print(
                        "Some processes failed excecution. Not all images may have processed. Error messages:"
                    )
                    for error in error_list:
                        if error is None:
                            continue
                        print(error)

        # Total Time
        logger.info("")
        logger.info("Total Time Taken = %f seconds" % (time.time() - st))


class Processor(object):
    def __init__(self, params, composite_tag=None, rank=0):
        self.params = params
        self.composite_tag = composite_tag

        # The convention is to put %s in the phil parameter to add a tag to
        # each output datafile. Save the initial templates here.
        self.experiments_filename_template = params.output.experiments_filename
        self.strong_filename_template = params.output.strong_filename
        self.indexed_filename_template = params.output.indexed_filename
        self.refined_experiments_filename_template = (
            params.output.refined_experiments_filename
        )
        self.integrated_filename_template = params.output.integrated_filename
        self.integrated_experiments_filename_template = (
            params.output.integrated_experiments_filename
        )

        debug_dir = os.path.join(params.output.output_dir, "debug")
        if not os.path.exists(debug_dir):
            try:
                os.makedirs(debug_dir)
            except OSError:
                pass  # due to multiprocessing, makedirs can sometimes fail
        assert os.path.exists(debug_dir)
        self.debug_file_path = os.path.join(debug_dir, "debug_%d.txt" % rank)
        write_newline = os.path.exists(self.debug_file_path)
        if write_newline:  # needed if the there was a crash
            self.debug_write("")

        if params.output.composite_output:
            assert composite_tag is not None
            from dxtbx.model.experiment_list import ExperimentList

            # self.all_strong_reflections = flex.reflection_table() # no composite strong pickles yet
            self.all_indexed_experiments = ExperimentList()
            self.all_indexed_reflections = flex.reflection_table()
            self.all_integrated_experiments = ExperimentList()
            self.all_integrated_reflections = flex.reflection_table()
            self.all_int_pickle_filenames = []
            self.all_int_pickles = []

            self.setup_filenames(composite_tag)

    def setup_filenames(self, tag):
        # before processing, set output paths according to the templates
        if (
            self.experiments_filename_template is not None
            and "%s" in self.experiments_filename_template
        ):
            self.params.output.experiments_filename = os.path.join(
                self.params.output.output_dir,
                self.experiments_filename_template % ("idx-" + tag),
            )
        if (
            self.strong_filename_template is not None
            and "%s" in self.strong_filename_template
        ):
            self.params.output.strong_filename = os.path.join(
                self.params.output.output_dir,
                self.strong_filename_template % ("idx-" + tag),
            )
        if (
            self.indexed_filename_template is not None
            and "%s" in self.indexed_filename_template
        ):
            self.params.output.indexed_filename = os.path.join(
                self.params.output.output_dir,
                self.indexed_filename_template % ("idx-" + tag),
            )
        if (
            self.refined_experiments_filename_template is not None
            and "%s" in self.refined_experiments_filename_template
        ):
            self.params.output.refined_experiments_filename = os.path.join(
                self.params.output.output_dir,
                self.refined_experiments_filename_template % ("idx-" + tag),
            )
        if (
            self.integrated_filename_template is not None
            and "%s" in self.integrated_filename_template
        ):
            self.params.output.integrated_filename = os.path.join(
                self.params.output.output_dir,
                self.integrated_filename_template % ("idx-" + tag),
            )
        if (
            self.integrated_experiments_filename_template is not None
            and "%s" in self.integrated_experiments_filename_template
        ):
            self.params.output.integrated_experiments_filename = os.path.join(
                self.params.output.output_dir,
                self.integrated_experiments_filename_template % ("idx-" + tag),
            )

    def debug_start(self, tag):
        import socket

        self.debug_str = "%s,%s" % (socket.gethostname(), tag)
        self.debug_str += ",%s,%s,%s\n"
        self.debug_write("start")

    def debug_write(self, string, state=None):
        from xfel.cxi.cspad_ana import cspad_tbx  # XXX move to common timestamp format

        ts = cspad_tbx.evt_timestamp()  # Now
        debug_file_handle = open(self.debug_file_path, "a")
        if string == "":
            debug_file_handle.write("\n")
        else:
            if state is None:
                state = "    "
            debug_file_handle.write(self.debug_str % (ts, state, string))
        debug_file_handle.close()

    def process_experiments(self, tag, experiments):

        if not self.params.output.composite_output:
            self.setup_filenames(tag)
        self.tag = tag
        self.debug_start(tag)

        if (
            not self.params.output.composite_output
            and self.params.output.experiments_filename
        ):

            experiments.as_json(self.params.output.experiments_filename)

        # Do the processing
        try:
            self.pre_process(experiments)
        except Exception as e:
            print("Error in pre-process", tag, str(e))
            self.debug_write("preprocess_exception", "fail")
            if not self.params.dispatch.squash_errors:
                raise
            return
        try:
            if self.params.dispatch.find_spots:
                self.debug_write("spotfind_start")
                observed = self.find_spots(experiments)
            else:
                print("Spot Finding turned off. Exiting")
                self.debug_write("data_loaded", "done")
                return
        except Exception as e:
            print("Error spotfinding", tag, str(e))
            self.debug_write("spotfinding_exception", "fail")
            if not self.params.dispatch.squash_errors:
                raise
            return
        try:
            if self.params.dispatch.index:
                if (
                    self.params.dispatch.hit_finder.enable
                    and len(observed)
                    < self.params.dispatch.hit_finder.minimum_number_of_reflections
                ):
                    print("Not enough spots to index", tag)
                    self.debug_write("not_enough_spots_%d" % len(observed), "stop")
                    return
                if (
                    self.params.dispatch.hit_finder.maximum_number_of_reflections
                    is not None
                ):
                    if (
                        self.params.dispatch.hit_finder.enable
                        and len(observed)
                        > self.params.dispatch.hit_finder.maximum_number_of_reflections
                    ):
                        print("Too many spots to index - Possibly junk", tag)
                        self.debug_write("too_many_spots_%d" % len(observed), "stop")
                        return
                self.debug_write("index_start")
                experiments, indexed = self.index(experiments, observed)
            else:
                print("Indexing turned off. Exiting")
                self.debug_write("spotfinding_ok_%d" % len(observed), "done")
                return
        except Exception as e:
            print("Couldn't index", tag, str(e))
            if not self.params.dispatch.squash_errors:
                raise
            self.debug_write("indexing_failed_%d" % len(observed), "stop")
            return
        self.debug_write("refine_start")
        try:
            experiments, indexed = self.refine(experiments, indexed)
        except Exception as e:
            print("Error refining", tag, str(e))
            self.debug_write("refine_failed_%d" % len(indexed), "fail")
            if not self.params.dispatch.squash_errors:
                raise
            return
        try:
            if self.params.dispatch.integrate:
                self.debug_write("integrate_start")
                integrated = self.integrate(experiments, indexed)
            else:
                print("Integration turned off. Exiting")
                self.debug_write("index_ok_%d" % len(indexed), "done")
                return
        except Exception as e:
            print("Error integrating", tag, str(e))
            self.debug_write("integrate_failed_%d" % len(indexed), "fail")
            if not self.params.dispatch.squash_errors:
                raise
            return
        self.debug_write("integrate_ok_%d" % len(integrated), "done")

    def pre_process(self, experiments):
        """ Add any pre-processing steps here """
        pass

    def find_spots(self, experiments):
        st = time.time()

        logger.info("*" * 80)
        logger.info("Finding Strong Spots")
        logger.info("*" * 80)

        # Find the strong spots
        observed = flex.reflection_table.from_observations(experiments, self.params)

        # Reset z coordinates for dials.image_viewer; see Issues #226 for details
        xyzobs = observed["xyzobs.px.value"]
        for i in range(len(xyzobs)):
            xyzobs[i] = (xyzobs[i][0], xyzobs[i][1], 0)
        bbox = observed["bbox"]
        for i in range(len(bbox)):
            bbox[i] = (bbox[i][0], bbox[i][1], bbox[i][2], bbox[i][3], 0, 1)

        if self.params.output.composite_output:
            pass  # no composite strong pickles yet
        else:
            # Save the reflections to file
            logger.info("\n" + "-" * 80)
            if self.params.output.strong_filename:
                self.save_reflections(observed, self.params.output.strong_filename)

        logger.info("")
        logger.info("Time Taken = %f seconds" % (time.time() - st))
        return observed

    def index(self, experiments, reflections):
        from dials.algorithms.indexing.indexer import Indexer

        st = time.time()

        logger.info("*" * 80)
        logger.info("Indexing Strong Spots")
        logger.info("*" * 80)

        params = copy.deepcopy(self.params)
        # don't do scan-varying refinement during indexing
        params.refinement.parameterisation.scan_varying = False

        if hasattr(self, "known_crystal_models"):
            known_crystal_models = self.known_crystal_models
        else:
            known_crystal_models = None

        if params.indexing.stills.method_list is None:
            idxr = Indexer.from_parameters(
                reflections,
                experiments,
                known_crystal_models=known_crystal_models,
                params=params,
            )
            idxr.index()
        else:
            indexing_error = None
            for method in params.indexing.stills.method_list:
                params.indexing.method = method
                try:
                    idxr = Indexer.from_parameters(
                        reflections, experiments, params=params
                    )
                    idxr.index()
                except Exception as e:
                    logger.info("Couldn't index using method %s" % method)
                    if indexing_error is None:
                        if e is None:
                            e = Exception("Couldn't index using method %s" % method)
                        indexing_error = e
                else:
                    indexing_error = None
                    break
            if indexing_error is not None:
                raise indexing_error

        indexed = idxr.refined_reflections
        experiments = idxr.refined_experiments

        if known_crystal_models is not None:

            filtered = flex.reflection_table()
            for idx in set(indexed["miller_index"]):
                sel = indexed["miller_index"] == idx
                if sel.count(True) == 1:
                    filtered.extend(indexed.select(sel))
            logger.info(
                "Filtered duplicate reflections, %d out of %d remaining"
                % (len(filtered), len(indexed))
            )
            print(
                "Filtered duplicate reflections, %d out of %d remaining"
                % (len(filtered), len(indexed))
            )
            indexed = filtered

        logger.info("")
        logger.info("Time Taken = %f seconds" % (time.time() - st))
        return experiments, indexed

    def refine(self, experiments, centroids):
        if self.params.dispatch.refine:
            from dials.algorithms.refinement import RefinerFactory

            st = time.time()

            logger.info("*" * 80)
            logger.info("Refining Model")
            logger.info("*" * 80)

            refiner = RefinerFactory.from_parameters_data_experiments(
                self.params, centroids, experiments
            )

            refiner.run()
            experiments = refiner.get_experiments()
            predicted = refiner.predict_for_indexed()
            centroids["xyzcal.mm"] = predicted["xyzcal.mm"]
            centroids["entering"] = predicted["entering"]
            centroids = centroids.select(refiner.selection_used_for_refinement())

            # Re-estimate mosaic estimates
            from dials.algorithms.indexing.nave_parameters import NaveParameters

            nv = NaveParameters(
                params=self.params,
                experiments=experiments,
                reflections=centroids,
                refinery=refiner,
                graph_verbose=False,
            )
            nv()
            acceptance_flags_nv = nv.nv_acceptance_flags
            centroids = centroids.select(acceptance_flags_nv)

        if self.params.output.composite_output:
            if (
                self.params.output.refined_experiments_filename
                or self.params.output.indexed_filename
            ):
                assert (
                    self.params.output.refined_experiments_filename is not None
                    and self.params.output.indexed_filename is not None
                )

                n = len(self.all_indexed_experiments)
                self.all_indexed_experiments.extend(experiments)
                for i, experiment in enumerate(experiments):
                    refls = centroids.select(centroids["id"] == i)
                    refls["id"] = flex.int(len(refls), n)
                    self.all_indexed_reflections.extend(refls)
                    n += 1
        else:
            # Dump experiments to disk
            if self.params.output.refined_experiments_filename:

                experiments.as_json(self.params.output.refined_experiments_filename)

            if self.params.output.indexed_filename:
                self.save_reflections(centroids, self.params.output.indexed_filename)

        if self.params.dispatch.refine:
            logger.info("")
            logger.info("Time Taken = %f seconds" % (time.time() - st))

        return experiments, centroids

    def integrate(self, experiments, indexed):
        st = time.time()

        logger.info("*" * 80)
        logger.info("Integrating Reflections")
        logger.info("*" * 80)

        indexed, _ = self.process_reference(indexed)

        # Get the integrator from the input parameters
        logger.info("Configuring integrator from input parameters")
        from dials.algorithms.profile_model.factory import ProfileModelFactory
        from dials.algorithms.integration.integrator import IntegratorFactory

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
            force_static=self.params.prediction.force_static,
        )
        predicted.match_with_reference(indexed)
        logger.info("")
        integrator = IntegratorFactory.create(self.params, experiments, predicted)

        # Integrate the reflections
        integrated = integrator.integrate()

        # correct integrated intensities for absorption correction, if necessary
        for abs_params in self.params.integration.absorption_correction:
            if abs_params.apply:
                if abs_params.algorithm == "fuller_kapton":
                    from dials.algorithms.integration.kapton_correction import (
                        multi_kapton_correction,
                    )
                elif abs_params.algorithm == "kapton_2019":
                    from dials.algorithms.integration.kapton_2019_correction import (
                        multi_kapton_correction,
                    )

                experiments, integrated = multi_kapton_correction(
                    experiments, integrated, abs_params.fuller_kapton, logger=logger
                )()

        if self.params.significance_filter.enable:
            from dials.algorithms.integration.stills_significance_filter import (
                SignificanceFilter,
            )
            from dxtbx.model.experiment_list import ExperimentList

            sig_filter = SignificanceFilter(self.params)
            filtered_refls = sig_filter(experiments, integrated)
            accepted_expts = ExperimentList()
            accepted_refls = flex.reflection_table()
            logger.info(
                "Removed %d reflections out of %d when applying significance filter"
                % (len(integrated) - len(filtered_refls), len(integrated))
            )
            for expt_id, expt in enumerate(experiments):
                refls = filtered_refls.select(filtered_refls["id"] == expt_id)
                if len(refls) > 0:
                    accepted_expts.append(expt)
                    refls["id"] = flex.int(len(refls), len(accepted_expts) - 1)
                    accepted_refls.extend(refls)
                else:
                    logger.info(
                        "Removed experiment %d which has no reflections left after applying significance filter"
                        % expt_id
                    )

            if len(accepted_refls) == 0:
                raise Sorry("No reflections left after applying significance filter")
            experiments = accepted_expts
            integrated = accepted_refls

        # Delete the shoeboxes used for intermediate calculations, if requested
        if self.params.integration.debug.delete_shoeboxes and "shoebox" in integrated:
            del integrated["shoebox"]

        if self.params.output.composite_output:
            if (
                self.params.output.integrated_experiments_filename
                or self.params.output.integrated_filename
            ):
                assert (
                    self.params.output.integrated_experiments_filename is not None
                    and self.params.output.integrated_filename is not None
                )

                n = len(self.all_integrated_experiments)
                self.all_integrated_experiments.extend(experiments)
                for i, experiment in enumerate(experiments):
                    refls = integrated.select(integrated["id"] == i)
                    refls["id"] = flex.int(len(refls), n)
                    self.all_integrated_reflections.extend(refls)
                    n += 1
        else:
            # Dump experiments to disk
            if self.params.output.integrated_experiments_filename:

                experiments.as_json(self.params.output.integrated_experiments_filename)

            if self.params.output.integrated_filename:
                # Save the reflections
                self.save_reflections(
                    integrated, self.params.output.integrated_filename
                )

        self.write_integration_pickles(integrated, experiments)
        from dials.algorithms.indexing.stills_indexer import (
            calc_2D_rmsd_and_displacements,
        )

        rmsd_indexed, _ = calc_2D_rmsd_and_displacements(indexed)
        log_str = "RMSD indexed (px): %f\n" % (rmsd_indexed)
        for i in range(6):
            bright_integrated = integrated.select(
                (
                    integrated["intensity.sum.value"]
                    / flex.sqrt(integrated["intensity.sum.variance"])
                )
                >= i
            )
            if len(bright_integrated) > 0:
                rmsd_integrated, _ = calc_2D_rmsd_and_displacements(bright_integrated)
            else:
                rmsd_integrated = 0
            log_str += (
                "N reflections integrated at I/sigI >= %d: % 4d, RMSD (px): %f\n"
                % (i, len(bright_integrated), rmsd_integrated)
            )

        for crystal_model in experiments.crystals():
            if hasattr(crystal_model, "get_domain_size_ang"):
                log_str += (
                    ". Final ML model: domain size angstroms: %f, half mosaicity degrees: %f"
                    % (
                        crystal_model.get_domain_size_ang(),
                        crystal_model.get_half_mosaicity_deg(),
                    )
                )

        logger.info(log_str)

        logger.info("")
        logger.info("Time Taken = %f seconds" % (time.time() - st))
        return integrated

    def write_integration_pickles(self, integrated, experiments, callback=None):
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
        if not hasattr(self.params.output, "integration_pickle"):
            return

        if self.params.output.integration_pickle is not None:
            from libtbx import easy_pickle
            from xfel.command_line.frame_extractor import ConstructFrame

            # Split everything into separate experiments for pickling
            for e_number, experiment in enumerate(experiments):
                e_selection = integrated["id"] == e_number
                reflections = integrated.select(e_selection)

                frame = ConstructFrame(reflections, experiment).make_frame()
                frame["pixel_size"] = experiment.detector[0].get_pixel_size()[0]

                if not hasattr(self, "tag") or self.tag is None:
                    try:
                        # if the data was a file on disc, get the path
                        event_timestamp = os.path.splitext(
                            experiments[0].imageset.paths()[0]
                        )[0]
                    except NotImplementedError:
                        # if the data is in memory only, check if the reader set a timestamp on the format object
                        event_timestamp = (
                            experiment.imageset.reader().get_format(0).timestamp
                        )
                    event_timestamp = os.path.basename(event_timestamp)
                    if event_timestamp.find("shot-") == 0:
                        event_timestamp = os.path.splitext(event_timestamp)[
                            0
                        ]  # micromanage the file name
                else:
                    event_timestamp = self.tag
                if hasattr(self.params.output, "output_dir"):
                    outfile = os.path.join(
                        self.params.output.output_dir,
                        self.params.output.integration_pickle
                        % (e_number, event_timestamp),
                    )
                else:
                    outfile = os.path.join(
                        os.path.dirname(self.params.output.integration_pickle),
                        self.params.output.integration_pickle
                        % (e_number, event_timestamp),
                    )

                if callback is not None:
                    callback(self.params, outfile, frame)

                if self.params.output.composite_output:
                    self.all_int_pickle_filenames.append(os.path.basename(outfile))
                    self.all_int_pickles.append(frame)
                else:
                    easy_pickle.dump(outfile, frame)

    def process_reference(self, reference):
        """ Load the reference spots. """
        if reference is None:
            return None, None
        st = time.time()
        assert "miller_index" in reference
        assert "id" in reference
        logger.info("Processing reference reflections")
        logger.info(" read %d strong spots" % len(reference))
        mask = reference.get_flags(reference.flags.indexed)
        rubbish = reference.select(mask == False)
        if mask.count(False) > 0:
            reference.del_selected(mask == False)
            logger.info(" removing %d unindexed reflections" % mask.count(True))
        if len(reference) == 0:
            raise Sorry(
                """
        Invalid input for reference reflections.
        Expected > %d indexed spots, got %d
      """
                % (0, len(reference))
            )
        mask = reference["miller_index"] == (0, 0, 0)
        if mask.count(True) > 0:
            rubbish.extend(reference.select(mask))
            reference.del_selected(mask)
            logger.info(" removing %d reflections with hkl (0,0,0)" % mask.count(True))
        mask = reference["id"] < 0
        if mask.count(True) > 0:
            raise Sorry(
                """
        Invalid input for reference reflections.
        %d reference spots have an invalid experiment id
      """
                % mask.count(True)
            )
        logger.info(" using %d indexed reflections" % len(reference))
        logger.info(" found %d junk reflections" % len(rubbish))
        logger.info(" time taken: %g" % (time.time() - st))
        return reference, rubbish

    def save_reflections(self, reflections, filename):
        """ Save the reflections to file. """
        st = time.time()
        logger.info("Saving %d reflections to %s" % (len(reflections), filename))
        reflections.as_pickle(filename)
        logger.info(" time taken: %g" % (time.time() - st))

    def finalize(self):
        """ Perform any final operations """
        if self.params.output.composite_output:
            if self.params.mp.composite_stride is not None:
                assert self.params.mp.method == "mpi"
                stride = self.params.mp.composite_stride

                from mpi4py import MPI

                comm = MPI.COMM_WORLD
                rank = comm.Get_rank()  # each process in MPI has a unique id, 0-indexed
                size = comm.Get_size()  # size: number of processes running in this job

                if rank % stride == 0:
                    subranks = [rank + i for i in range(1, stride) if rank + i < size]
                    for i in range(len(subranks)):
                        logger.info("Rank %d waiting for sender" % rank)
                        sender, indexed_experiments, indexed_reflections, integrated_experiments, integrated_reflections, int_pickles, int_pickle_filenames = comm.recv(
                            source=MPI.ANY_SOURCE
                        )
                        logger.info(
                            "Rank %d recieved data from rank %d" % (rank, sender)
                        )

                        if len(indexed_experiments) > 0:
                            indexed_reflections["id"] += len(
                                self.all_indexed_experiments
                            )
                            self.all_indexed_reflections.extend(indexed_reflections)
                            self.all_indexed_experiments.extend(indexed_experiments)

                        if len(integrated_experiments) > 0:
                            integrated_reflections["id"] += len(
                                self.all_integrated_experiments
                            )
                            self.all_integrated_reflections.extend(
                                integrated_reflections
                            )
                            self.all_integrated_experiments.extend(
                                integrated_experiments
                            )

                        self.all_int_pickles.extend(int_pickles)
                        self.all_int_pickle_filenames.extend(int_pickle_filenames)

                else:
                    destrank = (rank // stride) * stride
                    logger.info(
                        "Rank %d sending results to rank %d"
                        % (rank, (rank // stride) * stride)
                    )
                    comm.send(
                        (
                            rank,
                            self.all_indexed_experiments,
                            self.all_indexed_reflections,
                            self.all_integrated_experiments,
                            self.all_integrated_reflections,
                            self.all_int_pickles,
                            self.all_int_pickle_filenames,
                        ),
                        dest=destrank,
                    )

                    self.all_indexed_experiments = (
                        self.all_indexed_reflections
                    ) = (
                        self.all_integrated_experiments
                    ) = (
                        self.all_integrated_reflections
                    ) = self.all_int_pickles = self.all_integrated_reflections = []

            # Dump composite files to disk
            if (
                len(self.all_indexed_experiments) > 0
                and self.params.output.refined_experiments_filename
            ):

                self.all_indexed_experiments.as_json(
                    self.params.output.refined_experiments_filename
                )

            if (
                len(self.all_indexed_reflections) > 0
                and self.params.output.indexed_filename
            ):
                self.save_reflections(
                    self.all_indexed_reflections, self.params.output.indexed_filename
                )

            if (
                len(self.all_integrated_experiments) > 0
                and self.params.output.integrated_experiments_filename
            ):

                self.all_integrated_experiments.as_json(
                    self.params.output.integrated_experiments_filename
                )

            if (
                len(self.all_integrated_reflections) > 0
                and self.params.output.integrated_filename
            ):
                self.save_reflections(
                    self.all_integrated_reflections,
                    self.params.output.integrated_filename,
                )

            # Create a tar archive of the integration dictionary pickles
            if len(self.all_int_pickles) > 0 and self.params.output.integration_pickle:
                tar_template_integration_pickle = self.params.output.integration_pickle.replace(
                    "%d", "%s"
                )
                outfile = (
                    os.path.join(
                        self.params.output.output_dir,
                        tar_template_integration_pickle % ("x", self.composite_tag),
                    )
                    + ".tar"
                )
                tar = tarfile.TarFile(outfile, "w")
                for i, (fname, d) in enumerate(
                    zip(self.all_int_pickle_filenames, self.all_int_pickles)
                ):
                    string = StringIO(pickle.dumps(d, protocol=2))
                    info = tarfile.TarInfo(name=fname)
                    info.size = len(string.buf)
                    info.mtime = time.time()
                    tar.addfile(tarinfo=info, fileobj=string)
                tar.close()


if __name__ == "__main__":
    with dials.util.show_mail_on_error():
        script = Script()
        script.run()

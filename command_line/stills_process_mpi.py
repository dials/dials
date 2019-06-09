#!/usr/bin/env python
#
# LIBTBX_SET_DISPATCHER_NAME dials.stills_process_mpi

from __future__ import absolute_import, division, print_function

import copy
import glob
import logging
import os
from time import time

import libtbx.load_env
from libtbx.utils import Abort, Sorry
from dxtbx.model import Detector
from dials.util import log
from dials.command_line.dials_import import ManualGeometryUpdater
from dials.command_line.stills_process import Script as base_script
from dials.command_line.stills_process import do_import, phil_scope
from dials.command_line.stills_process import Processor
from dials.util.options import OptionParser

logger = logging.getLogger("dials.command_line.stills_process_mpi")

help_message = """
MPI derivative of dials.stills_process.  Only handle individual images, not HDF5
"""


class Script(base_script):
    """A class for running the script."""

    def __init__(self, comm):
        """MPI-aware constructor."""
        self.comm = comm
        self.rank = comm.Get_rank()  # each process in MPI has a unique id, 0-indexed
        self.size = comm.Get_size()  # size: number of processes running in this job

        # The script usage
        usage = (
            "usage: %s [options] [param.phil] mp.blob=<filepattern>"
            % libtbx.env.dispatcher_name
        )

        self.tag = None
        self.reference_detector = None

        # Create the parser
        self.parser = OptionParser(usage=usage, phil=phil_scope, epilog=help_message)

    def assign_work(self):

        """Execute the script."""

        if self.rank == 0:
            # Parse the command line
            params, options, all_paths = self.parser.parse_args(
                show_diff_phil=False, return_unhandled=True, quick_parse=True
            )

            # Check that all filenames have been entered as mp.blob
            assert all_paths == []
            assert params.mp.glob is not None
            # Log the diff phil
            diff_phil = self.parser.diff_phil.as_str()
            if diff_phil != "":
                logger.info("The following parameters have been modified:\n")
                logger.info(diff_phil)
                print(diff_phil)

            for item in params.mp.glob:
                all_paths += glob.glob(item)

            transmitted_info = dict(p=params, o=options, a=all_paths)
        else:
            transmitted_info = None

        transmitted_info = self.comm.bcast(transmitted_info, root=0)

        # Save the options
        self.options = transmitted_info["o"]
        self.params = transmitted_info["p"]
        all_paths = transmitted_info["a"]

        # Configure logging
        log.config(self.params.verbosity, info=None, debug=None)

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
        # Process the data
        assert self.params.mp.method == "mpi"

        basenames = [
            os.path.splitext(os.path.basename(filename))[0] for filename in all_paths
        ]
        tags = []
        for i, basename in enumerate(basenames):
            if basenames.count(basename) > 1:
                tags.append("%s_%05d" % (basename, i))
            else:
                tags.append(basename)
        iterable = zip(tags, all_paths)

        self.subset = [
            item for i, item in enumerate(iterable) if (i + self.rank) % self.size == 0
        ]
        print("DELEGATE %d of %d: %s" % (self.rank, self.size, self.subset[0:10]))

    def run(self):
        st = time()
        self.load_reference_geometry()

        update_geometry = ManualGeometryUpdater(self.params)

        # Import stuff
        # no preimport for MPI multifile specialization

        # Wrapper function
        def do_work(i, item_list):
            processor = Processor(copy.deepcopy(self.params), composite_tag="%04d" % i)
            for item in item_list:
                tag, filename = item

                experiments = do_import(filename)
                imagesets = experiments.imagesets()
                if len(imagesets) == 0 or len(imagesets[0]) == 0:
                    logger.info("Zero length imageset in file: %s" % filename)
                    return
                if len(imagesets) > 1:
                    raise Abort("Found more than one imageset in file: %s" % filename)
                if len(imagesets[0]) > 1:
                    raise Abort(
                        "Found a multi-image file. Run again with pre_import=True"
                    )

                if self.reference_detector is not None:
                    imagesets[0].set_detector(
                        Detector.from_dict(self.reference_detector.to_dict())
                    )

                update_geometry(imagesets[0])

                processor.process_experiments(tag, experiments)
            processor.finalize()

        # Process the data
        assert self.params.mp.method == "mpi"

        do_work(self.rank, self.subset)

        # Total Time
        logger.info("")
        logger.info("Total Time Taken = %f seconds" % (time() - st))


if __name__ == "__main__":
    from mpi4py import MPI

    comm = MPI.COMM_WORLD
    script = Script(comm)
    script.assign_work()
    comm.barrier()
    script.run()

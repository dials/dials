import dxtbx  # noqa: F401; import dependency to find HDF5 library

import os
import time
import numpy as np
import h5py
import glob
import bitshuffle.h5
import logging

from dxtbx.model.experiment_list import ExperimentListFactory
from dials.array_family import flex
from dials_algorithms_event_mode_ext import event_list

from dials_research.nexus.make_nxs import CBF2NexusWriter, CopyNexusStructure

logger = logging.getLogger("I2E")


def input_expt_names(expt_list):
    return sum(map(glob.glob, expt_list), [])


def get_input_image(imageset, i, msk=None):
    image = imageset.get_raw_data(i)[0]
    image.set_selected(~imageset.get_mask(i)[0], -1)
    if msk:
        image.set_selected(~msk, -2)
    return image


class event_generator:
    """
    A class to manage conversion of images to event data.
    """

    def __init__(self, image, num, exposure_time, pos_dset, time_dset):
        self._image = image
        self._num = num
        # Exposure time
        self._T = exposure_time
        # Output datasets
        self._d_position = pos_dset
        self._d_time = time_dset

    def write_to_file(self, events_position, events_time):
        n = events_time.size()

        P = self._d_position
        n_position = P.shape[0]
        P.resize(n_position + n, axis=0)
        P[-n:] = events_position.as_numpy_array()

        T = self._d_time
        n_time = T.shape[0]
        T.resize(n_time + n, axis=0)
        T[-n:] = events_time.as_numpy_array()

    def add_counts(self):
        tau0 = time.time()
        image = self._image
        sel = image > 0
        idx = sel.iselection()
        counts = image.select(idx)

        _pos, _t = event_list(self._num, self._T, idx, counts)
        num_counts = _pos.size()
        tau1 = time.time()
        logger.info(f"Total number of counts in this frame: {num_counts}")
        logger.info(f"Time taken for conversion: {tau1-tau0:0.2f}s")
        self.write_to_file(_pos, _t)


def images_to_events(params):
    t0 = time.time()

    # Define logger
    logdir = os.path.dirname(params.output.events)
    logging.basicConfig(
        filename=os.path.join(logdir, "images2events.log"),
        format="%(message)s",
        level="DEBUG",
    )

    # Get imageset from input
    expt = ExperimentListFactory.from_filenames(
        input_expt_names(params.input.experiments)
    )[0]
    imageset = expt.imageset

    # Get image range
    if params.input.image_range:
        n_img = params.input.image_range[-1]
    else:
        n_img = imageset.size()

    # Get "exposure time"
    exposure_time = params.input.exposure_time

    # Bad pixel mask
    # TODO TOFIX
    if params.input.mask:
        with h5py.File(params.input.mask[0], "r") as fh:
            msk = flex.bool(fh["bad_pixel_mask"][()])
    else:
        msk = None

    # Print out some information
    logger.info("Raw data: %s" % params.input.experiments)
    logger.info("Number of images to be converted: %d" % n_img)
    logger.info("Exposure time per frame: %d" % exposure_time)
    logger.info("")

    # Create output event data file
    logger.info("Events will be written to %s" % params.output.events)
    with h5py.File(params.output.events, "x") as fout:
        block_size = 0
        fout.create_dataset("time_per_frame", data=exposure_time)
        pos_dset = fout.create_dataset(
            "event_id",
            (0,),
            maxshape=(None,),
            dtype="i4",
            chunks=(100000,),
            compression=bitshuffle.h5.H5FILTER,
            compression_opts=(block_size, bitshuffle.h5.H5_COMPRESS_LZ4),
        )
        time_dset = fout.create_dataset(
            "event_time",
            (0,),
            maxshape=(None,),
            dtype="i4",
            chunks=(100000,),
            compression=bitshuffle.h5.H5FILTER,
            compression_opts=(block_size, bitshuffle.h5.H5_COMPRESS_LZ4),
        )

        # Generate events
        for j in range(n_img):
            logger.info("Converting image %d of %d" % (j, n_img))
            t1 = time.time()
            image = get_input_image(imageset, j, msk)
            event_generator(image, j, exposure_time, pos_dset, time_dset).add_counts()
            fout.flush()
            t2 = time.time()
            logger.info(f"Time taken: {t2 - t1:0.2f}s")
            logger.info("-" * 20)

    t3 = time.time()
    logger.info(f"Time taken for conversion: {t3 - t0:0.2f}s")
    # Get metadata from input file and write nxs
    logger.info("Writing experiment metadata to NeXus file...")
    # Check that the detector dimensions are saved correctly in input file
    # If not, prompt nexus writer to flip them
    fin, ext = os.path.splitext(params.input.experiments[0])
    if ext == ".h5" or ext == ".nxs":
        with h5py.File(params.input.experiments[0], "r") as fh:
            data_size = fh["entry/instrument/detector/module/data_size"][()]

        if np.all(data_size != imageset.get_raw_data(0)[0].all()):
            CopyNexusStructure(
                params.output.events,
                params.input.experiments[0],
                event_mode=True,
                flip=True,
            )
        else:
            CopyNexusStructure(
                params.output.events, params.input.experiments[0], event_mode=True,
            )
    else:
        CBF2NexusWriter(
            params.output.events,
            input_expt_names(params.input.experiments),
            event_mode=True,
        ).write_nexus_file()
    t4 = time.time()
    logger.info(f"Total processing time: {t4 - t0:0.2f}s")
    logger.info("=" * 20 + " EOF " + "=" * 20)


# TODO: save image_range information when it is applied
# Also save the corresponding scan axis values in the output

import dxtbx  # noqa: F401; import dependency to find HDF5 library

import os
import sys
import time
import numpy as np
import h5py
import bitshuffle.h5
import logging

from dials.array_family import flex
from dials_algorithms_event_mode_ext import image_coord

from make_nxs import CopyNexusStructure

CLOCK_FREQUENCY = int(6.4e8)
SHUTTER_OPEN = 0x840
SHUTTER_CLOSE = 0x880

logger = logging.getLogger("E2I")


def shutter_times(wdir):
    "Find timestamps for shutter open and close messages."
    "TO do this without calling the vds, I need to find which files that have the message. THose are not always the same but vary each experiment."
    "First thought, go through the directory and read all cue_id. This is not ideal once there are more modules"
    Topen = np.array([])
    Tclose = np.array([])
    for filename in os.listdir(wdir):
        ext = os.path.splitext(filename)[1]
        if ext == ".h5":
            if "vds" in filename:
                continue
            try:
                with h5py.File(os.path.join(wdir, filename), "r") as fh:
                    cues = fh["cue_id"][()]
                    if any(cues == SHUTTER_OPEN):
                        Topen = np.append(
                            Topen,
                            fh["cue_timestamp_zero"][
                                np.flatnonzero(cues == SHUTTER_OPEN)
                            ],
                        )
                    if any(cues == SHUTTER_CLOSE):
                        Tclose = np.append(
                            Tclose,
                            fh["cue_timestamp_zero"][
                                np.flatnonzero(cues == SHUTTER_CLOSE)
                            ],
                        )
            except KeyError:
                continue
    assert all(Topen) and all(Tclose), "Shutter signals do not match in the modules"
    "A better way to do this would probably be to look at the meta file and open only the files that have the first and last chunks of events. WIP."
    # for filename in os.listdir(wdir):
    #    if "meta" in filename:
    #        with h5py.File(filename, "r") as fh:
    #            for k in fh.keys():
    #                first = np.min(np.nonzero(fh[k][()]))
    #                last = np.max(np.nonzero(fh[k][()]))
    #                num = int(k.split("_")[-1]) + 1
    #                s = str(num).zfill(6)
    #                fname = meta.filename.replace("meta", s)
    # Open only the files with min and max idx and look for cues there.
    return Topen[0].astype(int), Tclose[0].astype(int)


def calc_num_frames(time_per_frame, t_i, t_f):
    " Returns the total number of images and the first and last frame number. "
    " time_per_frame is exposure_time(user input)* clock_freq"
    if (t_f - t_i) % time_per_frame == 0:
        n_frames = (t_f - t_i) // time_per_frame
    else:
        n_frames = (t_f - t_i) // time_per_frame + 1
    # first_frame = t_i // time_per_frame
    # last_frame = t_f // time_per_frame
    return n_frames


def decode_coordinates(event_pos, shape):
    "Returns decoded coordinates as flex.double()."
    "event_pos is passed as an array"
    "Timepix event_id dataset encodes the event coordinates as a 32bit word where the first 13 are x and the last 13 are y, with 0000 in between."
    x = event_pos & 0x1FFF
    y = event_pos >> 13
    xy = flex.double()
    for n in range(len(event_pos)):
        xy.append((y[n] + x[n] * shape[1]).astype(float))
    return xy


def histogram_image(coords, shape):
    " Returns the histogrammed data. Coords must be flex.double. "
    shape = tuple(map(int, shape))  # workaround. flex.grid does not work with np.int32
    dim = shape[0] * shape[1]
    H = flex.histogram(
        coords.as_double(), data_min=0.0, data_max=int(dim), n_slots=int(dim)
    )
    pixels = H.slots()
    pixels.reshape(flex.grid(shape))
    return pixels


class imager:
    """
    A class to convert pseudo-event mode data into images.
    """

    def __init__(self, expT, img_shape, output_dset):
        # Image shape
        self._shape = img_shape
        # Exposure time
        self._expT = int(expT)
        # Output
        self._dset = output_dset

    def write_to_file(self, image, n):
        start_dset = self._dset[n, :, :]
        self._dset[n, :, :] = start_dset + image.as_numpy_array()

    def _img_cache(self, n, pos, t):
        "Store events not yet assigned for next cycle"
        d = (t >= n * self._expT) & (t < (n + 1) * self._expT)
        idx = d.iselection()
        cp = pos.select(idx)
        ct = t.select(idx)
        return cp, ct

    def add_events(self, event_pos, event_time):
        chunk_size = event_time.chunks[0]
        num_events = event_time.len()
        if num_events % chunk_size == 0:
            chunk_count = num_events // chunk_size
        else:
            chunk_count = (num_events // chunk_size) + 1
        logger.info("%d chunks of %d size to be processed." % (chunk_count, chunk_size))

        cache_pos = flex.int()
        chache_time = flex.size_t()
        count = 0
        for n in range(chunk_count):
            tau0 = time.time()
            _pos = flex.int(event_pos[n * chunk_size : (n + 1) * chunk_size])
            _t = flex.size_t(
                event_time[n * chunk_size : (n + 1) * chunk_size].astype("<u8")
            )
            _pos = cache_pos.concatenate(_pos)
            _t = chache_time.concatenate(_t)

            start = flex.min(_t) // self._expT
            stop = (flex.max(_t) // self._expT) + 1
            for num in range(start, stop):
                if (num + 1 == stop) and (n + 1 != chunk_count):
                    cache_pps, cache_time = self._img_cache(num, _pos, _t)
                else:
                    coords = image_coord(num, self._expT, _pos, _t)
                    image = histogram_image(coords, self._shape)
                    self.write_to_file(image, num - 1)
            tau1 = time.time()
            logger.info(f"Time taken to process chunk {n}: {tau1-tau0:0.2f}s")
            events_written = _pos.size() - cache_pos.size()
            logger.info(f"Number of events written: {events_written}")
            count = count + events_written
        # Check that all events have been written
        if count != num_events:
            logger.warning(
                f"Not all events have been written! Number of missing events: {count - num_events}s"
            )


class timepix_imager:
    """
    A class to convert event data from timepix detector into images.
    """

    def __init__(self, t_i, t_f, expT, img_shape, output_dset):
        # Image shape
        self._shape = img_shape
        # Start and stop timestamps
        self._t0 = t_i
        self._t1 = t_f
        # Exposure time
        self._expT = expT
        # Output
        self._dset = output_dset

    def write_to_file(self, image, n):
        start_dset = self._dset[n, :, :]
        self._dset[n, :, :] = start_dset + image.as_numpy_array()

    def add_events(self, event_pos, event_time):
        chunk_size = event_time.chunks[0]
        num_events = event_time.len()
        if num_events % chunk_size == 0:
            chunk_count = num_events // chunk_size
        else:
            chunk_count = (num_events // chunk_size) + 1
        logger.info("%d chunks of %d size to be processed." % (chunk_count, chunk_size))

        for n in range(chunk_count):
            tau0 = time.time()
            _pos = flex.int(event_pos[n * chunk_size : (n + 1) * chunk_size])
            _t = flex.long(event_time[n * chunk_size : (n + 1) * chunk_size])
            # Decode coordinates
            _pos = decode_coordinates(_pos.as_numpy_array(), self._shape)
            # Discard all events recorded before and after shutter signals.
            indices = ((_t > self._t0) & (_t < self._t1)).iselection()
            # indices = (_t > self._t0).iselection()
            _pos = _pos.select(indices)
            _t = _t.select(indices)
            # Convert timestamp to frame number
            num = (_t - self._t0) / self._expT
            for _n in sorted(set(num)):
                ev = (num == _n).iselection()
                coords = _pos.select(ev)
                image = histogram_image(coords, self._shape)
                self.write_to_file(image, _n)
            tau1 = time.time()
            logger.info(f"Time taken to process chunk {n}: {tau1-tau0:0.2f}s")
            events_written = _pos.size() - _pos.count(0)
            logger.info(f"Number of events written: {events_written}")


def events_to_images(event_data: h5py.File, exposure_time, image_file: h5py.File):
    t0 = time.time()

    # Define logger
    logdir = os.path.dirname(image_file.filename)
    logging.basicConfig(
        filename=os.path.join(logdir, "events2images.log"),
        format="%(message)s",
        level="DEBUG",
    )

    # Get some information out of the Nexus file
    shape = event_data["entry/instrument/detector/module/data_size"][()]
    event_dir = os.path.dirname(event_data.filename)

    # Check if event data comes from Timepix detector
    try:
        description = event_data["entry/instrument/detector/description"][()].decode()
        logger.info("Detector description: %d" % description)
    except KeyError:
        logger.info("Detector description is missing from NeXus tree.")
        logger.info(
            "Looking for specific datasets to understand if data coming from Timepix detector..."
        )
        if "cue_id" in event_data["entry/data"].keys():
            description = "Timepix"
            logger.info("Assigned detector description: %s", description)
        else:
            description = "X"
            logger.info("No assigned detector description.")

    if description == "Timepix":
        # do something
        t_i, t_f = shutter_times(event_dir)
        time_per_frame = int(exposure_time * CLOCK_FREQUENCY)
        n_frames = calc_num_frames(time_per_frame, t_i, t_f)
    else:
        # do something else
        time_per_frame = event_data["entry/data/time_per_frame"][()]
        # Look for number of frames
        for k in event_data["entry/data"].keys():
            if "transformation_type" in event_data["entry/data"][k].attrs.keys():
                n_frames = event_data["entry/data"][k].shape[0]

    # Print out some information
    logger.info("Raw data directory: %s" % event_dir)
    logger.info("Detector dimensions (%d, %d)" % (shape[0], shape[1]))
    logger.info("Exposure time per frame: %d s" % exposure_time)
    logger.info("Number of images to be written: %d" % n_frames)
    logger.info("")

    # Create output image dataset
    logger.info("Images will be written to %s" % image_file.filename)
    logger.info("-" * 20)
    block_size = 0  # let Bitshuffle choose this
    dset = image_file.create_dataset(
        "data",
        shape=(n_frames, shape[0], shape[1]),
        dtype="i4",
        chunks=(1, shape[0], shape[1]),
        compression=bitshuffle.h5.H5FILTER,
        compression_opts=(block_size, bitshuffle.h5.H5_COMPRESS_LZ4),
    )

    t1 = time.time()
    if description == "Timepix":
        for filename in os.listdir(event_dir):
            ext = os.path.splitext(filename)[1]
            if ext == ".h5":
                if "meta" in filename or "vds" in filename:
                    continue
                try:
                    t2 = time.time()
                    with h5py.File(os.path.join(event_dir, filename), "r") as fh:
                        logger.info("Converting events from file %s" % filename)
                        event_id = fh["event_id"]
                        event_time = fh["event_time_offset"]
                        timepix_imager(
                            t_i, t_f, time_per_frame, shape, dset
                        ).add_events(event_id, event_time)
                        t3 = time.time()
                        logger.info(f"Time taken: {t3 - t2:0.2f}s")
                        logger.info("-" * 20)
                except OSError:  # for now, in case there's a subdirectory
                    continue
    else:
        logger.info("Converting events from file ... ")
        event_id = event_data["entry/data/event_id"]
        event_time = event_data["entry/data/event_time"]
        imager(time_per_frame, shape, dset).add_events(event_id, event_time)
        t4 = time.time()
        logger.info(f"Time taken: {t4 - t1:0.2f}s")
        logger.info("-" * 20)

    # Copy metadata from input file
    t5 = time.time()
    logger.info("")
    logger.info("Writing metadata to NeXus file...")
    CopyNexusStructure(image_file.filename, event_data.filename)
    t6 = time.time()
    logger.info(f"Time taken: {t6 - t5:0.2f}s")
    logger.info(f"Total processing time: {t6 - t0:0.2f}s")
    logger.info("=" * 20 + " EOF " + "=" * 20)


if __name__ == "__main__":
    event_data = sys.argv[1]
    exposure_time = sys.argv[2]
    image_file = sys.argv[3]
    with h5py.File(event_data, "r") as fin, h5py.File(image_file, "x") as fout:
        events_to_images(fin, exposure_time, fout)

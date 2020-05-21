from __future__ import division, print_function

from dxtbx.model.experiment_list import ExperimentListFactory
from dials.array_family import flex
import dials_algorithms_event_mode_ext

import dxtbx  # noqa: F401; import dependency to find HDF5 library

import h5py
import glob


class images_to_events(object):
    """
    Class to wrap/manage conversion of images to events.
    """

    def __init__(self, params):
        self._params = params
        self._expt = ExperimentListFactory.from_filenames(self.input_expt_names())[0]
        self._imageset = self._expt.imageset

        if self._params.input.image_range:
            self._n_img = self._params.input.image_range[-1]
        else:
            self._n_img = self._imageset.size()

        if self._params.input.mask:
            with h5py.File(self._params.input.mask[0], "r") as fh:
                self._msk = flex.bool(fh["bad_pixel_mask"][()])
        else:
            self._msk = None

        if self._params.output.events:
            self._fout = h5py.File(self._params.output.events, "x")
            self._d_position = self._fout.create_dataset(
                "position",
                (0,),
                maxshape=(None,),
                dtype="i4",
                chunks=(100000,),
                compression="lzf",
            )
            self._d_time = self._fout.create_dataset(
                "time",
                (0,),
                maxshape=(None,),
                dtype="i4",
                chunks=(100000,),
                compression="lzf",
            )
        else:
            print("No output file")

    def input_expt_names(self):
        return sum(map(glob.glob, self._params.input.experiments), [])

    def output(self, events_position, events_time):
        n = events_time.size()

        P = self._d_position
        n_position = P.shape[0]
        P.resize(n_position + n, axis=0)
        P[-n:] = events_position.as_numpy_array()

        T = self._d_time
        n_time = T.shape[0]
        T.resize(n_time + n, axis=0)
        T[-n:] = events_time.as_numpy_array()

        self._fout.flush()

    def img_input(self, i):
        image = self._imageset.get_raw_data(i)[0]
        image.set_selected(~self._imageset.get_mask(i)[0], -1)
        if self._msk is not None:
            image.set_selected(~self._msk, -2)
        return image

    def events(self, i):
        img = self.input_img(i)

        sel = img > 0
        idx = sel.iselection()
        counts = img.select(idx)

        pos, t = dials_algorithms_event_mode_ext.event_list(i, idx, counts)
        return pos, t

    def run(self):
        n = self._n_img
        for j in range(n):
            pos, t = self.events(j)
            if self._params.output.events:
                self.output(pos, t)
        if self._params.output.events:
            self._fout.close()

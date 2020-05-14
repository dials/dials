from __future__ import division, print_function

import dxtbx  # noqa: F401; import dependency to find HDF5 library

import sys
import h5py

from dials.array_family import flex
import dials_research_events


class images_to_events(object):
    """
    Class to wrap/manage conversion of images to events.
    """

    def __init__(self, h5_in, h5_out, h5_pix_mask, h5_meta):
        # h5_in -> nxs
        self._fin = h5py.File(h5_in, "r")
        self._data = self._fin["/entry/data/data"]
        self._shape = self._data.shape

        # h5_out
        self._fout = h5py.File(h5_out, "x")
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

        # pix_mask
        with h5py.File(h5_pix_mask, "r") as f1:
            self._pix_msk = flex.bool(f1["bad_pixel_mask"][()])

        # det_mask
        with h5py.File(h5_meta, "r") as f2:
            m = flex.int(f2["mask"][()])
            self._det_msk = m == 0

        # module shape, gap size (x,y)
        self._module = (1028, 512)
        self._gap = (12, 38)

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

    def img_input(self, n):
        image = flex.int(self._data[n].astype("i4"))
        image.set_selected(~self._det_msk, -1)
        image.set_selected(~self._pix_msk, -2)
        return image

    def events(self, n):
        image = self.img_input(n)

        mod_x = self._module[0]
        mod_y = self._module[1]
        gap_x = self._gap[0]
        gap_y = self._gap[1]

        # 8x4 m
        for i in range(8):
            for j in range(4):
                mod = image[
                    i * (mod_y + gap_y) : i * (mod_y + gap_y) + mod_y,
                    j * (mod_x + gap_x) : j * (mod_x + gap_x) + mod_x,
                ]
                sel = mod > 0
                idx = sel.iselection()
                counts = mod.select(idx)

                pos, t = dials_research_events.event_list(n, idx, counts)
                self.output(pos, t)

        return

    def run(self):
        for k in range(0, self._shape[0], 32):
            for _k in range(k, k + 32):
                if _k == self._shape[0]:
                    break

                self.events(_k)

        self._fout.close()


if __name__ == "__main__":
    i2e = images_to_events(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]).run()

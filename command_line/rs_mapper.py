#!/usr/bin/env python
#
# dials.rs_mapper.py
#
#  Copyright (C) 2014 Takanori Nakane
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import absolute_import, division, print_function

import math

import dials.algorithms.rs_mapper as recviewer
import dials.util
from cctbx import sgtbx, uctbx
from iotbx import ccp4_map, phil
from scitbx.array_family import flex

help_message = """
This program reconstructs reciprocal space from diffraction images. The orientation matrix is not necessary; only diffraction geometry is required.

This program is inteded to help detection and visualization of pathologies such as multiple-lattice, twinning, modulation, diffuse scattering and high background. It is also useful for education.

Examples::

  dials.rs_mapper image1.cbf

  dials.rs_mapper image_00*.cbf

  dials.rs_mapper imported.expt

"""

phil_scope = phil.parse(
    """
rs_mapper
  .short_caption = Reciprocal space mapper
{
  map_file = rs_mapper_output.ccp4
    .type = path
    .optional = False
    .multiple= False
    .short_caption = Map file
  max_resolution = 6
    .type = float
    .optional = True
    .short_caption = Resolution limit
  grid_size = 192
    .type = int
    .optional = True
  reverse_phi = False
    .type = bool
    .optional = True
}
""",
    process_includes=True,
)


class Script(object):
    def __init__(self):
        """Initialise the script."""
        from dials.util.options import OptionParser

        # The script usage
        usage = (
            "usage: dials.rs_mapper map_file=output.ccp4 [max_resolution=6] [grid_size=192] "
            "[reverse_phi=False] [param.phil] "
            "{image1.file [image2.file ...]} | imported.expt"
        )

        # Initialise the base class
        self.parser = OptionParser(
            usage=usage, phil=phil_scope, epilog=help_message, read_experiments=True
        )

    def run(self):
        from dials.util.options import flatten_experiments

        # Parse the command line
        params, options = self.parser.parse_args(show_diff_phil=True)

        if not params.rs_mapper.map_file:
            raise RuntimeError("Please specify output map file (map_file=)")
        else:
            self.map_file = params.rs_mapper.map_file

        # Ensure we have either a data block or an experiment list
        self.experiments = flatten_experiments(params.input.experiments)
        if len(self.experiments) != 1:
            self.parser.print_help()
            print("Please pass either an experiment list\n")
            return

        self.reverse_phi = params.rs_mapper.reverse_phi
        self.grid_size = params.rs_mapper.grid_size
        self.max_resolution = params.rs_mapper.max_resolution

        self.grid = flex.double(
            flex.grid(self.grid_size, self.grid_size, self.grid_size), 0
        )
        self.cnts = flex.int(
            flex.grid(self.grid_size, self.grid_size, self.grid_size), 0
        )

        for experiment in self.experiments:
            self.process_imageset(experiment.imageset)

        recviewer.normalize_voxels(self.grid, self.cnts)
        # Let's use 1/(100A) as the unit so that the absolute numbers in the
        # "cell dimensions" field of the ccp4 map are typical for normal
        # MX maps. The values in 1/A would give the "cell dimensions" around
        # or below 1 and some MX programs would not handle it well.
        box_size = 100 * 2.0 / self.max_resolution
        uc = uctbx.unit_cell((box_size, box_size, box_size, 90, 90, 90))
        ccp4_map.write_ccp4_map(
            self.map_file,
            uc,
            sgtbx.space_group("P1"),
            (0, 0, 0),
            self.grid.all(),
            self.grid,
            flex.std_string(["cctbx.miller.fft_map"]),
        )

    def process_imageset(self, imageset):
        rec_range = 1 / self.max_resolution

        panel = imageset.get_detector()[0]
        beam = imageset.get_beam()
        s0 = beam.get_s0()
        pixel_size = panel.get_pixel_size()
        xlim, ylim = imageset.get_raw_data(0)[0].all()

        # cache transformation
        xy = recviewer.get_target_pixels(panel, s0, xlim, ylim, self.max_resolution)

        s1 = panel.get_lab_coord(xy * pixel_size[0])  # FIXME: assumed square pixel
        s1 = s1 / s1.norms() * (1 / beam.get_wavelength())
        S = s1 - s0

        for i in range(len(imageset)):
            axis = imageset.get_goniometer().get_rotation_axis()
            osc_range = imageset.get_scan(i).get_oscillation_range()
            print("Oscillation range: %.1f - %.1f" % (osc_range[0], osc_range[1]))
            angle = (osc_range[0] + osc_range[1]) / 2 / 180 * math.pi
            if not self.reverse_phi:  # FIXME: ???
                angle *= -1
            rotated_S = S.rotate_around_origin(axis, angle)
            recviewer.fill_voxels(
                imageset.get_raw_data(i)[0],
                self.grid,
                self.cnts,
                rotated_S,
                xy,
                rec_range,
            )


if __name__ == "__main__":
    with dials.util.show_mail_on_error():
        script = Script()
        script.run()

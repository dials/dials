from __future__ import annotations

import math

from cctbx import sgtbx, uctbx
from iotbx import ccp4_map, phil
from scitbx.array_family import flex

import dials.algorithms.rs_mapper as recviewer
import dials.util
from dials.util import Sorry
from dials.util.options import ArgumentParser, flatten_experiments

help_message = """
This program reconstructs reciprocal space from diffraction images. The orientation matrix is not necessary; only diffraction geometry is required.

This program is intended to help detection and visualization of pathologies such as multiple-lattice, twinning, modulation, diffuse scattering and high background. It is also useful for education.

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
  map_file = rs_mapper_output.mrc
    .type = path
    .optional = False
    .multiple= False
    .short_caption = Output map file
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
  ignore_mask = False
    .type = bool
    .optional = True
    .short_caption = Ignore masks from dxtbx class
}
""",
    process_includes=True,
)


class Script:
    def __init__(self):
        """Initialise the script."""

        # The script usage
        usage = (
            "usage: dials.rs_mapper map_file=output.ccp4 [max_resolution=6] [grid_size=192] "
            "[reverse_phi=False] [param.phil] "
            "{image1.file [image2.file ...]} | imported.expt"
        )

        # Initialise the base class
        self.parser = ArgumentParser(
            usage=usage, phil=phil_scope, epilog=help_message, read_experiments=True
        )

    def run(self, args=None):
        # Parse the command line
        params, options = self.parser.parse_args(args, show_diff_phil=True)

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
        self.ignore_mask = params.rs_mapper.ignore_mask

        self.grid = flex.double(
            flex.grid(self.grid_size, self.grid_size, self.grid_size), 0
        )
        self.counts = flex.int(
            flex.grid(self.grid_size, self.grid_size, self.grid_size), 0
        )

        for experiment in self.experiments:
            self.process_imageset(experiment.imageset)

        recviewer.normalize_voxels(self.grid, self.counts)

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

        if len(imageset.get_detector()) != 1:
            raise Sorry("This program does not support multi-panel detectors.")

        panel = imageset.get_detector()[0]
        beam = imageset.get_beam()
        s0 = beam.get_s0()
        pixel_size = panel.get_pixel_size()
        xlim, ylim = imageset.get_raw_data(0)[0].all()
        if pixel_size[0] != pixel_size[1]:
            raise Sorry("This program does not support non-square pixels.")

        # cache transformation
        xy = recviewer.get_target_pixels(panel, s0, xlim, ylim, self.max_resolution)
        s1 = panel.get_lab_coord(xy * pixel_size[0])
        s1 = s1 / s1.norms() * (1 / beam.get_wavelength())
        S = s1 - s0

        for i in range(len(imageset)):
            axis = imageset.get_goniometer().get_rotation_axis()
            osc_range = imageset.get_scan(i).get_oscillation_range()
            print(f"Oscillation range: {osc_range[0]:.2f} - {osc_range[1]:.2f}")
            angle = (osc_range[0] + osc_range[1]) / 2 / 180 * math.pi
            if not self.reverse_phi:
                # the pixel is in S AFTER rotation. Thus we have to rotate BACK.
                angle *= -1
            rotated_S = S.rotate_around_origin(axis, angle)

            data = imageset.get_raw_data(i)[0]
            if not self.ignore_mask:
                mask = imageset.get_mask(i)[0]
                data.set_selected(~mask, 0)

            recviewer.fill_voxels(
                data, self.grid, self.counts, rotated_S, xy, rec_range
            )


@dials.util.show_mail_handle_errors()
def run(args=None):
    script = Script()
    script.run(args)


if __name__ == "__main__":
    run()

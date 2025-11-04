from __future__ import annotations

import concurrent.futures
import logging
import math

import numpy as np

import libtbx
from cctbx import sgtbx, uctbx
from iotbx import ccp4_map, phil
from scitbx.array_family import flex

import dials.algorithms.rs_mapper as recviewer
import dials.util
import dials.util.log
from dials.util import Sorry
from dials.util.options import ArgumentParser, flatten_experiments
from dials.util.system import CPU_COUNT

help_message = """
This program reconstructs reciprocal space from diffraction images. The orientation matrix is not necessary; only diffraction geometry is required.

This program is intended to help detection and visualization of pathologies such as multiple-lattice, twinning, modulation, diffuse scattering and high background. It is also useful for education.

Examples::

  dials.rs_mapper image1.cbf

  dials.rs_mapper image_00*.cbf

  dials.rs_mapper imported.expt
"""

# Define a logger
logger = logging.getLogger("dials.rs_mapper")

phil_scope = phil.parse(
    """
rs_mapper
  .short_caption = Reciprocal space mapper
{
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
  nproc = Auto
    .help = "Number of processes over which to split the calculation. If set to"
            "Auto, DIALS will choose automatically."
    .type = int(value_min=1)
    .expert_level = 1
}
output
{
  map_file = rs_mapper_output.mrc
    .type = path
    .optional = False
    .multiple= False
    .short_caption = Output map file

  log = "dials.rs_mapper.log"
    .type = str
}
""",
    process_includes=True,
)


def process_block(
    block, imageset, i_panel, grid_size, reverse_phi, S, ignore_mask, xy, rec_range
):
    grid = flex.double(flex.grid(grid_size, grid_size, grid_size), 0)
    counts = flex.int(flex.grid(grid_size, grid_size, grid_size), 0)

    axis = imageset.get_goniometer().get_rotation_axis()
    for i in block:
        osc_range = imageset.get_scan(i).get_oscillation_range()

        angle = (osc_range[0] + osc_range[1]) / 2 / 180 * math.pi
        if not reverse_phi:
            # the pixel is in S AFTER rotation. Thus we have to rotate BACK.
            angle *= -1
        rotated_S = S.rotate_around_origin(axis, angle)

        data = imageset.get_raw_data(i)[i_panel]
        if not ignore_mask:
            mask = imageset.get_mask(i)[i_panel]
            data.set_selected(~mask, 0)

        recviewer.fill_voxels(data, grid, counts, rotated_S, xy, rec_range)

    return grid, counts


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
            usage=usage,
            phil=phil_scope,
            epilog=help_message,
            read_experiments=True,
            read_experiments_from_images=True,
        )

    def run(self, args=None):
        # Parse the command line
        params, options = self.parser.parse_args(args, show_diff_phil=True)

        # Configure the logging.
        dials.util.log.config(options.verbose, logfile=params.output.log)

        if not params.output.map_file:
            raise RuntimeError("Please specify output map file (map_file=)")
        else:
            self.map_file = params.output.map_file

        self.experiments = flatten_experiments(params.input.experiments)

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

        self.nproc = params.rs_mapper.nproc
        if self.nproc is libtbx.Auto:
            self.nproc = CPU_COUNT
            logger.info(f"Setting nproc={self.nproc}")

        for i_expt, experiment in enumerate(self.experiments):
            logger.info(f"Calculation for experiment {i_expt}")
            for i_panel in range(len(experiment.detector)):
                grid, counts = self.process_imageset(experiment.imageset, i_panel)

                self.grid += grid
                self.counts += counts

        recviewer.normalize_voxels(self.grid, self.counts)

        # Let's use 1/(100A) as the unit so that the absolute numbers in the
        # "cell dimensions" field of the ccp4 map are typical for normal
        # MX maps. The values in 1/A would give the "cell dimensions" around
        # or below 1 and some MX programs would not handle it well.
        box_size = 100 * 2.0 / self.max_resolution
        uc = uctbx.unit_cell((box_size, box_size, box_size, 90, 90, 90))
        logger.info(f"Saving map to {self.map_file}")
        ccp4_map.write_ccp4_map(
            self.map_file,
            uc,
            sgtbx.space_group("P1"),
            (0, 0, 0),
            self.grid.all(),
            self.grid,
            flex.std_string(["cctbx.miller.fft_map"]),
        )

    def process_imageset(self, imageset, i_panel):
        rec_range = 1 / self.max_resolution

        beam = imageset.get_beam()
        s0 = beam.get_s0()

        panel = imageset.get_detector()[i_panel]
        pixel_size = panel.get_pixel_size()
        nfast, nslow = panel.get_image_size()

        if pixel_size[0] != pixel_size[1]:
            raise Sorry("This program does not support non-square pixels.")

        # cache transformation
        xy = recviewer.get_target_pixels(panel, s0, nfast, nslow, self.max_resolution)
        s1 = panel.get_lab_coord(xy * pixel_size[0])
        s1 = s1 / s1.norms() * (1 / beam.get_wavelength())
        S = s1 - s0

        # Split imageset into up to nproc blocks of at least 10 images
        nblocks = min(self.nproc, int(math.ceil(len(imageset) / 10)))
        blocks = np.array_split(range(len(imageset)), nblocks)
        blocks = [block.tolist() for block in blocks]

        logger.info(f"Calculation for panel {i_panel} split over {len(blocks)} blocks")
        header = ["Block", "Oscillation range (°)"]
        scan = imageset.get_scan()
        rows = [
            [
                f"{i + 1}",
                f"{scan.get_angle_from_array_index(block[0]):.2f} - {scan.get_angle_from_array_index(block[-1] + 1):.2f}",
            ]
            for i, block in enumerate(blocks)
        ]
        logger.info(dials.util.tabulate(rows, header, numalign="right") + "\n")

        if len(blocks) == 1:
            results = [
                process_block(
                    blocks[0],
                    imageset,
                    i_panel,
                    self.grid_size,
                    self.reverse_phi,
                    S,
                    self.ignore_mask,
                    xy,
                    rec_range,
                ),
            ]
        else:
            with concurrent.futures.ProcessPoolExecutor(
                max_workers=len(blocks)
            ) as pool:
                results = [
                    pool.submit(
                        process_block,
                        block,
                        imageset,
                        i_panel,
                        self.grid_size,
                        self.reverse_phi,
                        S,
                        self.ignore_mask,
                        xy,
                        rec_range,
                    )
                    for block in blocks
                ]
            results = [e.result() for e in results]

        grid, counts = results[0]
        for g, c in results[1:]:
            grid += g
            counts += c

        return grid, counts


@dials.util.show_mail_handle_errors()
def run(args=None):
    script = Script()
    script.run(args)


if __name__ == "__main__":
    run()

from __future__ import annotations

import copy
import json
import logging
import random
import sys

import numpy as np

import iotbx.phil
import libtbx
from dxtbx import flumpy
from scitbx.array_family import flex

import dials.util
from dials.algorithms.beam_position.compute_beam_position import compute_beam_position
from dials.algorithms.beam_position.helper_functions import (
    print_progress,
)
from dials.algorithms.beam_position.labelit_method import (
    discover_better_experimental_model,
)
from dials.util import Sorry, log
from dials.util.options import ArgumentParser, reflections_and_experiments_from_files
from dials.util.slice import slice_reflections
from dials.util.system import CPU_COUNT

logger = logging.getLogger("dials.command_line.search_beam_position")

help_message = """

A function to find beam center from diffraction images

The default method (based on the work of Sauter et al., J. Appl. Cryst.
37, 399-409 (2004)) is using the results from spot finding.

Example:
  dials.search_beam_position imported.expt strong.refl

Other methods are based on horizontal and vertical projection, and only
require an imported experiment.

Example:
  dials.search_beam_position method=midpoint imported.exp

More information about the projection methods can be found at
https://autoed.readthedocs.io/en/latest/pages/beam_position_methods.html

"""

phil_scope = iotbx.phil.parse(
    """

method = default midpoint maximum inversion

default {
    nproc = Auto
      .type = int(value_min=1)
    plot_search_scope = False
      .type = bool
    max_cell = None
      .type = float
      .help = "Known max cell (otherwise will compute from spot positions)"
    image_range = None
      .help = "The range of images to use in indexing. Number of arguments"
        "must be a factor of two. Specifying \"0 0\" will use all images"
        "by default. The given range follows C conventions"
        "(e.g. j0 <= j < j1)."
      .type = ints(size=2)
      .multiple = True
    max_reflections = 10000
      .type = int(value_min=1)
      .help = "Maximum number of reflections to use in the search for better"
              "experimental model. If the number of input reflections is "
              "greater then a random subset of reflections will be used."
    mm_search_scope = 4.0
      .help = "Global radius of origin offset search."
      .type = float(value_min=0)
    wide_search_binning = 2
      .help = "Modify the coarseness of the wide grid search for "
              "the beam centre."
      .type = float(value_min=0)
    n_macro_cycles = 1
      .type = int
      .help = "Number of macro cycles for an iterative beam centre search."
    d_min = None
      .type = float(value_min=0)
    seed = 42
      .type = int(value_min=0)
}

projection {

    method_x = midpoint maximum inversion
    .type = str
    .help = "The projection method along the x-axis."

    method_y = midpoint maximum inversion
    .type = str
    .help = "The projection method along the y-axis."

    plot = True
    .type = bool
    .help = "Plot the diffraction image with the computed beam center."

    bar = True
    .type = bool
    .help = "Print progress bar."

    exclude_pixel_range_x = None
    .type = ints
    .multiple = True
    .help = "List of comma-separated pairs of numbers specifying pixel ranges "
            "in the x direction to exclude from projection to the y-axis "
            "(e.g., exclude_pixel_range_x=20,350,700,800 would exclude ranges "
            "20-350 and 700-800). Indexing assumes Python (or C) "
            "conventions. The first pixel has an index 0, the last pixel has "
            "an index N-1 (here N is the number of pixels along the x-axis). "
            "The last pixel in the range is not included, e.g., '0,N' would "
            "exclude the entire range, while '0,3' would exclude pixels 0, 1, "
            "and 2."

    exclude_pixel_range_y = None
    .type = ints
    .multiple = True
    .help = "List of pixel ranges to exclude from projection to the y-axis. "
            "See `exclude_pixel_range_x` for more details."

    per_image = False
    .type = bool
    .help = "Compute the beam position for each image. Otherwise, compute the "
            "beam position for a single (average) image."

    color_cutoff = None
    .type = float
    .help = "The maximum of the colorbar range in the plotted beam position "
            "figure. Use this option to adjust the visibility of the plotted "
            "diffraction image."

    midpoint {

        exclude_intensity_percent = 0.01
        .type = float
        .help = "Order all pixels by intensity and discard this percentage"
                "from the top (by setting them to zero)."

        intersection_range = (0.3, 0.9, 0.01)
        .type = floats
        .help = "Compute midpoints in this range (start, end, step)."

        convolution_width = 80
        .type = int
        .help = "Convolution kernel width used for smoothing (in pixels)."

        dead_pixel_range_x = None
        .type = ints
        .multiple = True
        .help = "List of comma-separated pairs of numbers specifying pixel "
                "ranges in the x direction to exclude from midpoint "
                "calculation. Indexing assumes the same rules as with "
                "the exclude_pixel_range_x."

        dead_pixel_range_y = None
        .type = ints
        .multiple = True
        .help = "List of comma-separated pairs of numbers specifying pixel "
                "ranges in the y direction to exclude from midpoint "
                "calculation. Indexing assumes the same rules as with "
                "exclude_pixel_range_y."

        intersection_min_width = 10
        .type = int
        .help = "Do not consider midpoint intersections below this width."

    }

    maximum {

        bad_pixel_threshold = None
        .type = int
        .help = "Set all pixels above this value to zero."

        n_convolutions = 1
        .type = int
        .help = "The number of smoothing convolutions."

        convolution_width = 1
        .type = int
        .help = "Convolution kernel width used for smoothing (in pixels)."

        bin_width = 20
        .type = int
        .help = "The width of the averaging bin used to find the region "
                "of max intensity (pixels)."

        bin_step = 10
        .type = int
        .help = "Distance in pixels between neighboring bins used to find the "
                "region of maximal intensity."

    }

    inversion {

        bad_pixel_threshold = None
        .type = int
        .help = "Set all pixels above this value to zero."

        guess_position = None
        .type = ints(size=2)
        .help = "Initial guess for the beam position (x, y) in pixels. "
                "If not supplied, it will be set to the center of the image."

        inversion_window_width = 400
        .type = int
        .help = "Do profile inversion within this window (in pixels)."

        background_cutoff = None
        .type = int
        .help = "Set all the pixels with intensity above this value to zero."

        convolution_width = 1
        .type = int
        .help = "Convolution kernel width used for smoothing (in pixels)."
    }

}

output {
  experiments = optimised.expt
    .type = path
  log = "dials.search_beam_position.log"
    .type = str
  json = "beam_positions.json"
    .type = str
}
"""
)


@dials.util.show_mail_handle_errors()
def run(args=None):
    usage = """
    dials.search_beam_position imported.expt strong.refl
    dials.search_beam_position method=midpoint imported.exp
    dials.search_beam_position method_x=midpoint method_y=maximum imported.exp
    """

    try:
        parser = ArgumentParser(
            usage=usage,
            phil=phil_scope,
            read_experiments=True,
            read_reflections=True,
            check_format=True,
            epilog=help_message,
        )

        params, options = parser.parse_args(args, show_diff_phil=False)
    except Sorry:
        parser = ArgumentParser(
            usage=usage,
            phil=phil_scope,
            read_experiments=True,
            read_reflections=True,
            check_format=False,
            epilog=help_message,
        )

        params, options = parser.parse_args(args, show_diff_phil=False)

    reflections, experiments = reflections_and_experiments_from_files(
        params.input.reflections, params.input.experiments
    )

    cond_01 = len(params.projection.method_x) == 1
    cond_02 = len(params.projection.method_y) == 1

    # Configure the logging
    log.config(logfile=params.output.log)

    if params.method[0] == "default" and not (cond_01 or cond_02):
        if len(experiments) == 0 or len(reflections) == 0:
            parser.print_help()
            exit(0)

        # Log the diff phil
        diff_phil = parser.diff_phil.as_str()
        if diff_phil != "":
            logger.info("The following parameters have been modified:\n")
            logger.info(diff_phil)

        if params.default.seed is not None:
            flex.set_random_seed(params.default.seed)
            random.seed(params.default.seed)

        if params.default.nproc is libtbx.Auto:
            params.default.nproc = CPU_COUNT

        imagesets = experiments.imagesets()
        # Split all the refln tables by ID, corresponding to
        # the respective imagesets
        reflections = [
            refl_unique_id
            for refl in reflections
            for refl_unique_id in refl.split_by_experiment_id()
        ]

        assert len(imagesets) > 0
        # assert len(reflections) == len(imagesets)
        if len(reflections) != len(imagesets):
            log.info("No reflections found.")
            sys.exit()

        if (
            params.default.image_range is not None
            and len(params.default.image_range) > 0
        ):
            reflections = [
                slice_reflections(refl, params.default.image_range)
                for refl in reflections
            ]

        for i in range(params.default.n_macro_cycles):
            if params.default.n_macro_cycles > 1:
                logger.info("Starting macro cycle %i", i + 1)
            experiments = discover_better_experimental_model(
                experiments,
                reflections,
                params.default,
                nproc=params.default.nproc,
                d_min=params.default.d_min,
                mm_search_scope=params.default.mm_search_scope,
                wide_search_binning=params.default.wide_search_binning,
                plot_search_scope=params.default.plot_search_scope,
            )
            logger.info("")

        msg = "Saving optimised experiments to %s" % params.output.experiments
        logger.info(msg)
        experiments.as_file(params.output.experiments)

    else:  # Other methods (midpoint, maximum, inversion)
        new_experiments = copy.deepcopy(experiments)

        if len(new_experiments) == 0:
            parser.print_help()
            sys.exit(0)

        # Compute the starting positions for all imagesets
        imageset_positions = {}
        start_index = 0
        for exp in new_experiments:
            imageset_positions[exp.imageset] = start_index
            start_index += exp.imageset.size()

        json_output = []

        # Map imagesets to detectors
        detector_to_imagesets = {d: [] for d in new_experiments.detectors()}
        for expt in new_experiments:
            detector_to_imagesets[expt.detector].append(expt.imageset)

        for det_index, detector in enumerate(detector_to_imagesets.keys()):
            average_image = None
            image_counter = 0
            positions_x = []
            positions_y = []

            n_imagesets = len(detector_to_imagesets[detector])

            for jj, imageset in enumerate(detector_to_imagesets[detector]):
                start_index = imageset_positions[imageset]
                end_index = start_index + imageset.size()

                for i, ii in enumerate(range(start_index, end_index)):
                    if params.projection.bar and not params.projection.per_image:

                        print_progress(
                            image_index=i,
                            n_images=abs(end_index - start_index),
                            set_index=jj,
                            n_sets=n_imagesets,
                        )

                    index = np.uint64(i)
                    image = imageset.get_corrected_data(index)
                    image = flumpy.to_numpy(image[0])

                    mask = imageset.get_mask(0)
                    mask = flumpy.to_numpy(mask[0])

                    image[mask == 0] = 0

                    if params.projection.per_image:
                        x, y = compute_beam_position(
                            image, params, image_counter, det_index
                        )

                        msg = f"[Detector {int(det_index)}] "
                        msg += f"image {int(image_counter):05d}: "
                        msg += f"{x:.2f} {y:.2f}"
                        logger.info(msg)

                        if x is not None and y is not None:
                            positions_x.append(x)
                            positions_y.append(y)
                            json_output.append(
                                (
                                    int(det_index),
                                    int(image_counter),
                                    float(x),
                                    float(y),
                                )
                            )

                        else:
                            log.warn("WARNING: beam position not found.")

                    else:
                        if average_image is None:
                            average_image = image
                        else:
                            average_image = average_image + image

                    image_counter += 1

            if params.projection.per_image:
                x = np.array(positions_x).mean()
                y = np.array(positions_y).mean()
                msg = 40 * "=" + "\n"
                msg += f"[Detector {det_index}] "
                msg += "Average (per image) beam position: \n"
                msg += f" {x:.2f} {y:.2f}"
                logger.info(msg)
            else:
                average_image = average_image / image_counter
                x, y = compute_beam_position(image, params, detector_index=det_index)
                msg = f"[Detector {det_index}] "
                msg += "Beam position from the average image: \n "
                msg += f" {x:.2f} {y:.2f}"
                logger.info(msg)

            if x is not None and y is not None:
                for panel in detector:
                    xy = panel.pixel_to_millimeter((float(x), float(y)))

                    beam = experiments[0].beam
                    beam_direction = np.array(beam.get_s0())

                    norm = np.linalg.norm(beam_direction)
                    beam_direction = beam_direction / norm

                    distance = panel.get_distance()
                    beam_vec = beam_direction * distance

                    fast = np.array(panel.get_fast_axis())
                    slow = np.array(panel.get_slow_axis())

                    new_origin = beam_vec - fast * xy[0] - slow * xy[1]
                    panel.set_frame(fast, slow, new_origin)

        msg = "Saving optimised experiments to %s" % params.output.experiments
        logger.info(msg)
        new_experiments.as_file(params.output.experiments)

        if params.projection.per_image:
            with open(params.output.json, "w") as json_file:
                json.dump(json_output, json_file, indent=4)


if __name__ == "__main__":
    run()

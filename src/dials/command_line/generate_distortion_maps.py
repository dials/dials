from __future__ import annotations

import logging
import math
import pickle

from iotbx import phil
from scitbx import matrix

from dials.algorithms.image.distortion import CreateEllipticalDistortionMaps
from dials.array_family import flex
from dials.util import Sorry, log, show_mail_handle_errors
from dials.util.options import ArgumentParser, flatten_experiments
from dials.util.version import dials_version

logger = logging.getLogger("dials.command_line.generate_distortion_maps")

help_message = """

Generate dx.pickle, dy.pickle distortion maps for a detector model picked up
from an image file or experiment.expt. These maps can be used to
represent distortion within the millimetre to pixel mapping

Examples::

  dials.generate_distortion_maps image_001.cbf dx=0.5 dy=1.5
  dials.generate_distortion_maps models.expt mode=ellipse phi=0 l2=0.95
"""

scope = phil.parse(
    """
  mode = *translate ellipse
    .type = choice

  translate
    .help = "Options for applying separate translation offsets to each panel"
  {
    dx = 0.0
      .type = floats
      .help = "Translation in pixels to be applied along the fast direction."
              "The number of values supplied should equal the number of panels."
    dy = 0.0
      .type = floats
      .help = "Translation in pixels to be applied along the slow direction."
              "The number of values supplied should equal the number of panels."
  }

  ellipse
    .help = "Options for correcting for elliptical distortion of images."
            "Defaults set for correction of datasets published in"
            "https://doi.org/10.1107/S2059798317010348"
  {
    phi = -21.0
      .type = float
      .help = "Acute angle of one principal axis of the ellipse from the fast"
              "axis of the first panel of the detector"
    l1 = 1.0
      .type = float
      .help = "Scale factor for first axis of the ellipse"

    l2 = 0.956
      .type = float
      .help = "Scale factor for second axis of the ellipse"

    centre_xy = 33.2475,33.2475
      .type = floats(size=2)
      .help = "Centre of the ellipse in millimetres along fast, slow of the"
              "first panel"
  }

  output {
    x_map = dx.pickle
      .type = str
    y_map = dy.pickle
      .type = str
    log = dials.generate_distortion_maps.log
      .type = str
  }
"""
)


def make_dx_dy_translate(imageset, dx, dy):
    images = imageset.indices()
    image = imageset[images[0]]

    if (len(dx) != len(image)) or (len(dx) != len(image)):
        raise Sorry(
            "Please provide separate translations for each panel of the detector"
        )

    distortion_map_x = []
    distortion_map_y = []

    for block, shift_x, shift_y in zip(image, dx, dy):
        distortion_map_x.append(flex.double(flex.grid(block.focus()), shift_x))
        distortion_map_y.append(flex.double(flex.grid(block.focus()), shift_y))

    return tuple(distortion_map_x), tuple(distortion_map_y)


def ellipse_matrix_form(phi, l1, l2):
    """Return the matrix for the quadratic form describing the oblique ellipse
    where the first axis makes an angle phi with the X axis and the scale factors
    for the axes are l1 and l2.
    See https://www.le.ac.uk/users/dsgp1/COURSES/TOPICS/quadrat.pdf"""
    deg2rad = math.pi / 180.0
    phi *= deg2rad
    cphi = math.cos(phi)
    sphi = math.sin(phi)

    a11 = l1 * cphi**2 + l2 * sphi**2
    a12 = (l2 - l1) * sphi * cphi
    a21 = a12
    a22 = l1 * sphi**2 + l2 * cphi**2

    assert a11 * a22 - 2 * a12 > 0.0

    return matrix.sqr((a11, a12, a21, a22))


def make_dx_dy_ellipse(imageset, phi, l1, l2, centre_xy):
    detector = imageset.get_detector()

    # Get matrix describing the elliptical distortion
    M = ellipse_matrix_form(phi, l1, l2)

    distortion_map_x = []
    distortion_map_y = []

    for panel in detector:
        map_maker = CreateEllipticalDistortionMaps(panel, M, centre_xy)
        distortion_map_x.append(map_maker.get_dx())
        distortion_map_y.append(map_maker.get_dy())

    distortion_map_x = tuple(distortion_map_x)
    distortion_map_y = tuple(distortion_map_y)

    return distortion_map_x, distortion_map_y


@show_mail_handle_errors()
def run(args=None):
    usage = "dials.generate_distortion_maps [options] image_*.cbf"

    parser = ArgumentParser(
        usage=usage,
        phil=scope,
        read_experiments=True,
        read_experiments_from_images=True,
        check_format=False,
        epilog=help_message,
    )

    params, options = parser.parse_args(args)

    # Configure the logging
    log.config(verbosity=options.verbose, logfile=params.output.log)

    logger.info(dials_version())

    # Log the diff phil
    diff_phil = parser.diff_phil.as_str()
    if diff_phil != "":
        logger.info("The following parameters have been modified:\n")
        logger.info(diff_phil)

    experiments = flatten_experiments(params.input.experiments)

    if len(experiments) == 0:
        parser.print_help()
        exit()

    assert len(experiments) == 1

    imagesets = experiments.imagesets()

    assert len(imagesets) == 1

    imageset = imagesets[0]

    if params.mode == "translate":
        op = params.translate
        logger.info(f"Generating translation map with dx={op.dx}, dy={op.dy}")
        dx, dy = make_dx_dy_translate(imageset, op.dx, op.dy)
    elif params.mode == "ellipse":
        op = params.ellipse
        logger.info(
            "Generating elliptical map with phi={}, l1={}, "
            "l2={}, centre_xy={},{}".format(op.phi, op.l1, op.l2, *op.centre_xy)
        )
        dx, dy = make_dx_dy_ellipse(imageset, op.phi, op.l1, op.l2, op.centre_xy)
    else:
        raise Sorry("Unrecognised mode")

    logger.info(f"Saving X distortion map to {params.output.x_map}")
    with open(params.output.x_map, "wb") as f:
        pickle.dump(dx, f, pickle.HIGHEST_PROTOCOL)

    logger.info(f"Saving Y distortion map to {params.output.y_map}")
    with open(params.output.y_map, "wb") as f:
        pickle.dump(dy, f, pickle.HIGHEST_PROTOCOL)


if __name__ == "__main__":
    run()

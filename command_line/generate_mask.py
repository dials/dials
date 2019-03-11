#!/usr/bin/env python
# -*- coding: utf-8 -*-

# generate_mask.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

"""
Mask images to remove unwanted pixels.

This program is used to generate mask to specify which pixels should be
considered "invalid" during spot finding and integration. It provides a few
options to create simple masks using the detector trusted range, or from
simple shapes or by setting different resolution ranges.

Examples::

  dials.generate_mask experiments.json border=5

  dials.generate_mask experiments.json \\
    untrusted.rectangle=50,100,50,100 \\
    untrusted.circle=200,200,100

  dials.generate_mask experiments.json resolution.d_max=2.00

"""

from __future__ import absolute_import, division, print_function

import sys
import logging
import six.moves.cPickle as pickle
from typing import Tuple, Optional

import dials.util
import dials.util.log
from scitbx.array_family import flex
from iotbx.phil import parse
from libtbx.phil import scope, scope_extract
import libtbx.load_env
from dials.util.options import OptionParser, flatten_experiments
from dials.util.masking import MaskGenerator
from dxtbx.format.image import ImageBool
from dxtbx.model.experiment_list import Experiment, ExperimentList, ExperimentListDumper

Masks = Tuple[flex.bool, ...]

help_message = __doc__

log = logging.getLogger('dials.generate_mask')

phil_scope = parse(
    """
    output {
        mask = mask.pickle
            .type = path
            .help = "Name of output mask file."
        experiments = None
            .type = path
            .help = "Name of output experiment list file.  If this is set, a copy of "
                    "the experiments, modified with the generated pixel masks, " 
                    "will be saved to this location."
        log = 'dials.generate_masks.log'
            .type = str
            .help = "The log filename."
    }

    include scope dials.util.masking.phil_scope
    
    verbosity = 1
        .type = int(value_min=0)
        .help = "The verbosity level."
    """,
    process_includes=True,
)


def generate_mask(experiments, params):
    # type: (ExperimentList, scope_extract) -> Tuple[Masks, Optional[ExperimentList]]
    """
    Generate a pixel mask for each image in an experiment.

    Use the masking parameters :param:`params` and an experiment in the experiment
    list :param:`experiments` to define pixel masks for the associated imagesets.
    The masks are generated using :mod:`dials.util.masking`.

    The masks will be saved to disk at the location specified by
    :attr:`params.output.mask`.  Optionally, if a path
    :attr:`params.output.experiments` is set, a modified copy of :param:`experiments`
    with the masks applied will be saved to that location.

    Args:
        experiments: An experiment list containing only one experiment.
        params: Masking parameters, having the structure defined in :data:`phil_scope`.

    Returns:
        A tuple containing the generated pixel masks.
        A copy of :param:`experiments` with the masks applied
            (optional, only returned if :attr:`params.output.experiments` is set).
    """
    imagesets = experiments.imagesets()
    if len(imagesets) > 1:
        # Check beams (for resolution) and detectors are equivalent in each case
        # otherwise the mask may not be appropriate across all imagesets
        detectors = experiments.detectors()
        beams = experiments.beams()
        for d in detectors[1:]:
            if not d.is_similar_to(detectors[0]):
                sys.exit(
                    "Multiple imagesets are present, but their detector models differ."
                )
        for b in beams[1:]:
            if not b.is_similar_to(beams[0]):
                sys.exit(
                    "Multiple imagesets are present, but their beam models differ."
                )

    imageset = imagesets[0]

    # Generate the mask
    generator = MaskGenerator(params)
    mask = generator.generate(imageset)

    # Save the mask to file
    log.info("Writing mask to %s", params.output.mask)
    with open(params.output.mask, "wb") as fh:
        pickle.dump(mask, fh)

    # Save the experiment list
    if params.output.experiments:
        for imageset in imagesets:
            imageset.external_lookup.mask.data = ImageBool(mask)
            imageset.external_lookup.mask.filename = params.output.mask

        log.info("Saving experiments to %s", params.output.experiments)
        dump = ExperimentListDumper(experiments)
        dump.as_file(params.output.experiments)

        return mask, experiments
    else:
        return mask, None


def run(phil=phil_scope, args=None):
    # type: (scope, list) -> None
    """
    Parse command-line arguments, run the script.

    Use the DIALS option parser to extract an experiment list and parameters.  Pass
    these to :func:`script`.  If :param:`args` is `None` (default), the option parser
    defaults to :data:`sys.argv[1:]`.

    Args:
        phil: PHIL scope for option parser.
        args: Arguments to parse.  Defaults to :data:`sys.argv[1:]`.
    """
    # Create the parser
    usage = "usage: %s [options] experiments.json" % libtbx.env.dispatcher_name
    parser = OptionParser(
        usage=usage, phil=phil, epilog=help_message, read_experiments=True
    )

    # Parse the command line arguments
    params, options = parser.parse_args(args=args, show_diff_phil=True)
    experiments = flatten_experiments(params.input.experiments)

    # Configure logging
    dials.util.log.config(params.verbosity, info=params.output.log)

    # Check number of args
    if len(experiments) == 0:
        parser.print_help()
        sys.exit(1)

    if len(experiments) != 1:
        sys.exit("Exactly 1 experiment must be specified.")

    # Run the script
    generate_mask(experiments, params)


if __name__ == "__main__":
    with dials.util.show_mail_on_error():
        run()

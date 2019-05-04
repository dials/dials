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

  dials.generate_mask datablock.json border=5

  dials.generate_mask datablock.json \\
    untrusted.rectangle=50,100,50,100 \\
    untrusted.circle=200,200,100

  dials.generate_mask datablock.json resolution.d_max=2.00

"""

from __future__ import absolute_import, division, print_function

import logging
import os.path
import sys

import six.moves.cPickle as pickle

import dials.util
import dials.util.log
import libtbx.load_env
import libtbx.phil as phil
from dials.util.masking import MaskGenerator
from dials.util.options import OptionParser, flatten_datablocks
from dxtbx.format.image import ImageBool
from dxtbx.datablock import DataBlock, DataBlockDumper
from scitbx.array_family import flex

try:
    from typing import List, Optional, Tuple

    Masks = List[Tuple[flex.bool, ...]]
except ImportError:
    pass

log = logging.getLogger("dials.generate_mask")

phil_scope = phil.parse(
    """
    output {
        mask = mask.pickle
            .type = path
            .help = "Name of output mask file."
        datablock = None
            .type = path
            .help = "Name of output datablock file.  If this is set, a copy of "
                    "the datablock, modified with the generated pixel masks, "
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


def generate_mask(datablock, params):
    # type: (DataBlock, phil.scope_extract) -> Tuple[Masks, Optional[DataBlock]]
    """
    Generate a pixel mask for each imageset in a datablock.

    Use the masking parameters :param:`params` to define a pixel mask for each of the
    imagesets in a datablock :param:`datablock`.  The masks are generated using
    :mod:`dials.util.masking`.

    The masks will be saved to disk at the location specified by
    :attr:`params.output.mask`.  If the datablock contains more than one imageset,
    multiple mask files will be produced, with filenames differentiated by an
    appended number.  Optionally, if a path :attr:`params.output.experiments` is
    set, a modified copy of :param:`datablock` with the masks applied will be saved
    to that location.

    Args:
        datablock: A datablock containing only one imageset.
        params: Masking parameters, having the structure defined in
            :data:`phil_scope`.

    Returns:
        A list of masks, one for each imageset.

        A copy of :param:`datablock` with the masks applied (optional,
        only returned if :attr:`params.output.datablock` is set).
    """
    imagesets = datablock.extract_imagesets()
    masks = []

    # Create output mask filenames
    num_imagesets = len(imagesets)
    if num_imagesets == 1:
        filenames = [params.output.mask]
    else:
        # If there is more than one imageset, append a number to each output filename
        name, ext = os.path.splitext(params.output.mask)
        pad = len(str(num_imagesets))
        filenames = [
            "{name}_{num:0{pad}}{ext}".format(name=name, num=i + 1, pad=pad, ext=ext)
            for i in range(num_imagesets)
        ]

    # Generate the mask
    generator = MaskGenerator(params)

    for imageset, filename in zip(imagesets, filenames):
        mask = generator.generate(imageset)
        masks.append(mask)

        # Save the mask to file
        log.info("Writing mask to %s", filename)
        with open(filename, "wb") as fh:
            pickle.dump(mask, fh)

        if params.output.datablock:
            # Apply the mask to the imageset
            imageset.external_lookup.mask.data = ImageBool(mask)
            imageset.external_lookup.mask.filename = filename

    if params.output.datablock:
        # Save the datablock
        log.info("Saving datablock to %s", params.output.datablock)
        dump = DataBlockDumper(datablock)
        dump.as_file(params.output.datablock)
    else:
        datablock = None

    return masks, datablock

def run(phil=phil_scope, args=None):
    # type: (phil.scope, List[str]) -> None
    """
    Parse command-line arguments, run the script.

    Uses the DIALS option parser to extract a datablock and parameters, then passes
    these to :func:`generate_mask`.

    Args:
        phil: PHIL scope for option parser.
        args: Arguments to parse. If None, :data:`sys.argv[1:]` will be used.
    """
    # Create the parser
    usage = "usage: %s [options] datablock.json" % libtbx.env.dispatcher_name
    parser = OptionParser(usage=usage, phil=phil, epilog=__doc__, read_datablocks=True)

    # Parse the command line arguments
    params, options = parser.parse_args(args=args, show_diff_phil=True)
    datablocks = flatten_datablocks(params.input.datablock)

    # Configure logging
    dials.util.log.config(params.verbosity, info=params.output.log)

    # Check number of args
    if len(datablocks) == 0:
        parser.print_help()
        sys.exit(1)
    if len(datablocks) > 1:
        log.warning('Please supply just one datablock file for masking.')
        sys.exit(1)

    # Run the script
    generate_mask(datablocks[0], params)


if __name__ == "__main__":
    with dials.util.show_mail_on_error():
        run()

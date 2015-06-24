#!/usr/bin/env python
#
# aaply_mask.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division
from iotbx.phil import parse

phil_scope = parse("""

  input {
    mask = None
      .type = str
      .help = "The mask filename"
  }

  output {
    datablock = datablock_with_mask.json
      .type = str
      .help = "Name of output datablock file"
  }
""", process_includes=True)


class Script(object):
  ''' A class to encapsulate the script. '''

  def __init__(self):
    ''' Initialise the script. '''
    from dials.util.options import OptionParser
    import libtbx.load_env

    # Create the parser
    usage = "usage: %s [options] datablock.json" % libtbx.env.dispatcher_name
    self.parser = OptionParser(
      usage=usage,
      phil=phil_scope,
      read_datablocks=True)

  def run(self):
    ''' Run the script. '''
    from dials.util.options import flatten_datablocks
    from dxtbx.datablock import DataBlockDumper
    from libtbx.utils import Abort
    import cPickle as pickle

    # Parse the command line arguments
    params, options = self.parser.parse_args(show_diff_phil=True)
    datablocks = flatten_datablocks(params.input.datablock)

    # Check number of args
    if len(datablocks) == 0:
      self.parser.print_help()
      return

    # Check the mask file is given
    if params.input.mask is None:
      self.parser.print_help()
      return

    # Check nbumber of datablocks
    if len(datablocks) != 1:
      raise Abort('exactly 1 datablock must be specified')

    # Get the imageset
    datablock = datablocks[0]
    imagesets = datablock.extract_imagesets()
    if len(imagesets) != 1:
      raise Abort('datablock must contain exactly 1 imageset')
    imageset = imagesets[0]

    # Set the lookup
    imageset.external_lookup.mask.filename = params.input.mask

    # Dump the datablock
    print "Writing datablock to %s" % params.output.datablock
    dump = DataBlockDumper(datablock)
    dump.as_json(filename=params.output.datablock)


if __name__ == '__main__':
  from dials.util import halraiser
  try:
    script = Script()
    script.run()
  except Exception as e:
    halraiser(e)

#!/usr/bin/env python
#
# generate_mask.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import absolute_import, division
from iotbx.phil import parse

help_message = '''

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

'''

phil_scope = parse("""
  output {
    mask = mask.pickle
      .type = str
      .help = "Name of output mask file"
  }

  include scope dials.util.masking.phil_scope
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
      epilog=help_message,
      read_datablocks=True)

  def run(self):
    ''' Run the script. '''
    from dials.util.masking import MaskGenerator
    from dials.util.options import flatten_datablocks
    from libtbx.utils import Sorry
    import cPickle as pickle
    from dials.util import log

    # Parse the command line arguments
    params, options = self.parser.parse_args(show_diff_phil=True)
    datablocks = flatten_datablocks(params.input.datablock)

    # COnfigu logging
    log.config()

    # Check number of args
    if len(datablocks) == 0:
      self.parser.print_help()
      return

    if len(datablocks) != 1:
      raise Sorry('exactly 1 datablock must be specified')
    datablock = datablocks[0]
    imagesets = datablock.extract_imagesets()
    if len(imagesets) != 1:
      raise Sorry('datablock must contain exactly 1 imageset')
    imageset = imagesets[0]

    # Generate the mask
    generator = MaskGenerator(params)
    mask = generator.generate(imageset)

    # Save the mask to file
    print "Writing mask to %s" % params.output.mask
    pickle.dump(mask, open(params.output.mask, "w"))

if __name__ == '__main__':
  from dials.util import halraiser
  try:
    script = Script()
    script.run()
  except Exception as e:
    halraiser(e)

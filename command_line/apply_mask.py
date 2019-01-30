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

from __future__ import absolute_import, division, print_function

from iotbx.phil import parse

help_message = '''

This program augments a datablock or experiment list file with a mask specified
by the user.

Examples::

  dials.apply_mask datablock.json input.mask=mask.pickle

  dials.apply_mask expriments.json input.mask=mask.pickle

'''

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

    experiments = experiments_with_mask.json
      .type = str
      .help = "Name of output experiments file"
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
      epilog=help_message,
      phil=phil_scope,
      read_experiments=True,
      read_datablocks=True)

  def run(self):
    ''' Run the script. '''
    from dials.util.options import flatten_datablocks
    from dials.util.options import flatten_experiments
    from dxtbx.datablock import DataBlockDumper
    from dxtbx.model.experiment_list import ExperimentListDumper
    from dials.util import Sorry

    # Parse the command line arguments
    params, options = self.parser.parse_args(show_diff_phil=True)

    # Check the mask file is given
    if params.input.mask is None:
      self.parser.print_help()
      return

    experiments = flatten_experiments(params.input.experiments)
    datablocks = flatten_datablocks(params.input.datablock)
    do_experiments = len(experiments) > 0
    if do_experiments and datablocks:
      self.parser.print_help()
      raise Sorry("Either a datablock or an experiment list may be provided"
                  " but not both together.")

    if datablocks:
      if len(datablocks) != 1:
        raise Sorry('exactly 1 datablock must be specified')
      datablock = datablocks[0]
      imagesets = datablock.extract_imagesets()
    elif do_experiments:
      imagesets = experiments.imagesets()
    else:
      raise Sorry("Either a datablock or an experiment list may be provided")

    # Get the imageset
    if len(imagesets) != 1:
      raise Sorry('A mask can be applied only to a single imageset')
    imageset = imagesets[0]

    # Set the lookup
    imageset.external_lookup.mask.filename = params.input.mask

    # Dump the datablock
    if datablocks:
      print("Writing datablock to %s" % params.output.datablock)
      dump = DataBlockDumper(datablock)
      dump.as_json(filename=params.output.datablock)
    else:
      print("Writing experiments to %s" % params.output.experiments)
      dump = ExperimentListDumper(experiments)
      dump.as_json(filename=params.output.experiments)


if __name__ == '__main__':
  from dials.util import halraiser
  try:
    script = Script()
    script.run()
  except Exception as e:
    halraiser(e)

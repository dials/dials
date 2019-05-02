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

import cPickle as pickle

from dxtbx.format.image import ImageBool
from iotbx.phil import parse

help_message = """

This program augments a datablock or experiement list .json file with one or more masks specified by the
user.  Its only function is to input the mask file paths to the .json file,
but means that the user does not have to edit the file by hand.

Crucially, the mask files must be provided in the same order as their corresponding
imagesets (sweeps) appear in the .json file.

Examples::

    dials.apply_mask datablock.json input.mask=mask.pickle

    dials.apply_mask experiments.json input.mask=mask1.pickle input.mask=mask2.pickle

"""

phil_scope = parse(
    """
        input {
            mask = None
                .multiple = True
                .type = str
                .help = "The mask filenames, one mask per imageset"
        }

  output {
    datablock = masked_datablock.json
      .type = str
      .help = "Name of output datablock file"

    experiments = masked_experiments.json
      .type = str
      .help = "Name of output experiments file"
  }
""",
    process_includes=True,
)


class Script(object):
    """ A class to encapsulate the script. """

    def __init__(self):
        """ Initialise the script. """
        from dials.util.options import OptionParser
        import libtbx.load_env

        # Create the parser
        usage = (
            "usage: %s datablock_or_experiments.json input.mask=mask.pickle"
            % libtbx.env.dispatcher_name
        )
        self.parser = OptionParser(
            usage=usage,
            epilog=help_message,
            phil=phil_scope,
            read_experiments=True,
            read_datablocks=True,
        )

    def run(self):
        """ Run the script. """
        from dials.util.options import flatten_datablocks
        from dials.util.options import flatten_experiments
        from dxtbx.datablock import DataBlockDumper
        from dxtbx.model.experiment_list import ExperimentListDumper
        from dials.util import Sorry

        # Parse the command line arguments
        params, options = self.parser.parse_args(show_diff_phil=True)
        experiments = flatten_experiments(params.input.experiments)
        datablocks = flatten_datablocks(params.input.datablock)

        # Check that a valid JSON file and at least one mask file have been provided
        if not ((experiments or datablocks) and params.input.mask):
            self.parser.print_help()
            return

        if experiments and datablocks:
            self.parser.print_help()
            raise Sorry(
                "Either a datablock or an experiment list may be provided "
                "but not both together."
            )

        if datablocks:
            if len(datablocks) != 1:
                raise Sorry("exactly 1 datablock must be specified")
            # Check number of datablocks matches the number of masks
            n_dblocks = len(datablocks)
            n_masks = len(params.input.mask)
            if n_dblocks != n_masks:
                raise Sorry(
                    "The number of masks provided must match the number of imagesets "
                    "(sweeps).\n"
                    "You have provided a datablock containing {} imageset(s).\n"
                    "You have provided {} mask file(s).".format(n_dblocks, n_masks)
                )
            datablock = datablocks[0]
            imagesets = datablock.extract_imagesets()
        elif experiments:
            # Check number of experiments matches the number of masks
            n_expts = len(experiments)
            n_masks = len(params.input.mask)
            if n_expts != n_masks:
                raise Sorry(
                    "The number of masks provided must match the number of imagesets "
                    "(sweeps).\n"
                    "You have provided an experiment list containing {} imageset(s).\n"
                    "You have provided {} mask file(s).".format(n_expts, n_masks)
                )
            imagesets = experiments.imagesets()
        else:
            raise Sorry("Either a datablock or an experiment list must be provided")

        for i, imageset in enumerate(imagesets):
            # Set the lookup
            with open(params.input.mask[i]) as f:
                mask = pickle.load(f)
            imageset.external_lookup.mask.filename = params.input.mask[i]
            imageset.external_lookup.mask.data = ImageBool(mask)

        # Dump the datablock or experiment list
        if datablocks:
            print("Writing datablock to %s" % params.output.datablock)
            dump = DataBlockDumper(datablock)
            dump.as_json(filename=params.output.datablock)
        else:
            print("Writing experiments to %s" % params.output.experiments)
            dump = ExperimentListDumper(experiments)
            dump.as_json(filename=params.output.experiments)


if __name__ == "__main__":
    from dials.util import halraiser

    try:
        script = Script()
        script.run()
    except Exception as e:
        halraiser(e)

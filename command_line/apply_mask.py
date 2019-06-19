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

from six.moves import cPickle as pickle

from dxtbx.format.image import ImageBool
from iotbx.phil import parse

help_message = """

This program augments a experiments JSON file with one or more masks specified by the
user.  Its only function is to input the mask file paths to the experiments JSON file,
but means that the user does not have to edit the experiments file by hand.

Crucially, the mask files must be provided in the same order as their corresponding
imagesets (sweeps) appear in the experiments JSON file.

Examples::

    dials.apply_mask models.expt input.mask=pixels.mask

    dials.apply_mask models.expt input.mask=pixels1.mask input.mask=pixels2.mask

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
            experiments = masked.expt
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
            "usage: %s models.expt input.mask=pixels.mask" % libtbx.env.dispatcher_name
        )
        self.parser = OptionParser(
            usage=usage, epilog=help_message, phil=phil_scope, read_experiments=True
        )

    def run(self):
        """ Run the script. """
        from dials.util.options import flatten_experiments
        from dxtbx.model.experiment_list import ExperimentListDumper
        from dials.util import Sorry

        # Parse the command line arguments
        params, options = self.parser.parse_args(show_diff_phil=True)
        experiments = flatten_experiments(params.input.experiments)

        # Check that an experiment list and at least one mask file have been provided
        if not (experiments and params.input.mask):
            self.parser.print_help()
            return

        # Check number of experiments
        n_expts = len(experiments)
        n_masks = len(params.input.mask)
        if n_expts != n_masks:
            raise Sorry(
                "The number of masks provided must match the number of imagesets "
                "(sweeps).\n"
                "You have provided an experiment list containing {} imageset(s).\n"
                "You have provided {} mask file(s).".format(n_expts, n_masks)
            )

        # Get the imageset
        imagesets = experiments.imagesets()

        for i, imageset in enumerate(imagesets):
            # Set the lookup
            with open(params.input.mask[i]) as f:
                mask = pickle.load(f)
            imageset.external_lookup.mask.filename = params.input.mask[i]
            imageset.external_lookup.mask.data = ImageBool(mask)

        # Dump the experiments
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

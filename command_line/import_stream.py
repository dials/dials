#!/usr/bin/env python
#
# import_stream.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is` distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import absolute_import, division, print_function

import logging

from libtbx.phil import parse
from dials.util import Sorry
from dxtbx.model.experiment_list import ExperimentListDumper
from dxtbx.model.experiment_list import ExperimentListFactory

logger = logging.getLogger("dials.command_line.import_stream")

help_message = """


"""

# Create the phil parameters

phil_scope = parse(
    """

  output {

    experiments = imported.expt
      .type = str
      .help = "The output JSON or pickle file"

    log = 'dials.import_stream.log'
      .type = str
      .help = "The log filename"

    debug_log = 'dials.import_stream.debug.log'
      .type = str
      .help = "The debug log filename"

    compact = False
      .type = bool
      .help = "For JSON output use compact representation"

    directory = auto
      .type = str
      .help = "The output directory for streaming data"

    image_template = "%05d.image"
      .type = str
      .help = "The image template"

  }

  verbosity = 0
    .type = int(value_min=0)
    .help = "The verbosity level"

  input {

    host = localhost
      .type = str
      .help = "The input host"

    port = 9999
      .type = int
      .help = "The input port"

  }

"""
)


class Script(object):
    """ Class to parse the command line options. """

    def __init__(self):
        """ Set the expected options. """
        from dials.util.options import OptionParser
        import libtbx.load_env

        # Create the option parser
        usage = "usage: %s [options]" % libtbx.env.dispatcher_name
        self.parser = OptionParser(
            usage=usage, sort_options=True, phil=phil_scope, epilog=help_message
        )

    def run(self):
        """ Parse the options. """
        from dials.util import log
        import libtbx
        from uuid import uuid4
        from dials.util.stream import ZMQStream, Decoder
        from os.path import join, exists
        import os
        import json

        # Parse the command line arguments in two passes to set up logging early
        params, options = self.parser.parse_args(show_diff_phil=False, quick_parse=True)

        # Configure logging
        log.config(
            params.verbosity, info=params.output.log, debug=params.output.debug_log
        )
        from dials.util.version import dials_version

        logger.info(dials_version())

        # Parse the command line arguments completely
        params, options = self.parser.parse_args(show_diff_phil=False)

        # Log the diff phil
        diff_phil = self.parser.diff_phil.as_str()
        if diff_phil != "":
            logger.info("The following parameters have been modified:\n")
            logger.info(diff_phil)

        # Check a stream is given
        if params.input.host is None:
            raise Sorry("An input host needs to be given")

        # Check the directory
        if params.output.directory is None:
            raise Sorry("An output directory needs to be given")
        elif params.output.directory is libtbx.Auto:
            params.output.directory = "/dev/shm/dials-%s" % uuid4()

        # Make the output directory
        if exists(params.output.directory):
            raise Sorry('Directory "%s" already exists' % (params.output.directory))

        # Make the directory
        os.mkdir(params.output.directory)

        # Initialise the stream
        stream = ZMQStream(params.input.host, params.input.port)
        decoder = Decoder(params.output.directory, params.output.image_template)
        imageset = None
        while True:

            # Get the frames from zmq
            frames = stream.receive()

            # Decode the frames
            obj = decoder.decode(frames)

            # Process the object
            if obj.is_header():
                filename = join(params.output.directory, "metadata.json")
                with open(filename, "w") as outfile:
                    json.dump(obj.header, outfile)
                imageset = obj.as_imageset(filename)
                experiments = ExperimentListFactory.from_imageset(imageset)
                self.write_experiments(experiments, params)
            elif obj.is_image():
                assert imageset is not None
                filename = join(
                    params.output.directory, params.output.image_template % obj.count
                )
                with open(filename, "wb") as outfile:
                    outfile.write(obj.data)
                filename = join(
                    params.output.directory,
                    "%s.info" % (params.output.image_template % obj.count),
                )
                with open(filename, "w") as outfile:
                    json.dump(obj.info, outfile)
            elif obj.is_endofseries():
                assert imageset is not None
                break
            else:
                raise RuntimeError("Unknown object")

        # Close the stream
        stream.close()

    def write_experiments(self, experiments, params):
        """
        Output the experiments to file.

        """
        if params.output.experiments:
            logger.info("-" * 80)
            logger.info("Writing experiments to %s" % params.output.experiments)
            dump = ExperimentListDumper(experiments)
            dump.as_file(params.output.experiments, compact=params.output.compact)


if __name__ == "__main__":
    from dials.util import halraiser

    try:
        script = Script()
        script.run()
    except Exception as e:
        halraiser(e)

#!/usr/bin/env python
#
# dials.find_shared_models.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#

from __future__ import absolute_import, division, print_function

import logging

logger = logging.getLogger("dials.command_line.find_shared_models")

help_message = """

This program attempts to find sets of images with shared models

Examples::

  dials.find_shared_models models.expt

  dials.find_shared_models experiments1.expt experiments2.expt

"""

# Set the phil scope
from libtbx.phil import parse

phil_scope = parse(
    """
  output {
    log = 'dials.find_shared_models.log'
      .type = str
      .help = "The log filename"

    debug_log = 'dials.find_shared_models.debug.log'
      .type = str
      .help = "The debug log filename"
  }

  verbosity = 0
    .type = int
    .help = "The verbosity level"
""",
    process_includes=True,
)


class Script(object):
    """A class for running the script."""

    def __init__(self):
        """Initialise the script."""
        from dials.util.options import OptionParser
        import libtbx.load_env

        # The script usage
        usage = (
            "usage: %s [options] [param.phil] "
            "models.expt " % libtbx.env.dispatcher_name
        )

        # Initialise the base class
        self.parser = OptionParser(
            usage=usage, phil=phil_scope, epilog=help_message, read_experiments=True
        )

    def run(self):
        """Execute the script."""
        from dials.util.options import flatten_experiments
        from dxtbx.imageset import ImageSweep
        from time import time
        from dials.util import log
        import datetime

        start_time = time()

        # Parse the command line
        params, options = self.parser.parse_args(show_diff_phil=False)

        # Configure the logging
        log.config(
            params.verbosity, info=params.output.log, debug=params.output.debug_log
        )

        from dials.util.version import dials_version

        logger.info(dials_version())

        # Log the diff phil
        diff_phil = self.parser.diff_phil.as_str()
        if diff_phil != "":
            logger.info("The following parameters have been modified:\n")
            logger.info(diff_phil)

        # Ensure we have a data block
        experiments = flatten_experiments(params.input.experiments)
        if len(experiments) == 0:
            self.parser.print_help()
            return

        # Get the list of sweeps
        sweeps = []
        for experiment in experiments:
            if isinstance(experiment.imageset, ImageSweep):
                sweeps.append(experiment.imageset)
        logger.info("Number of sweeps = %d" % len(sweeps))

        # Sort the sweeps by timestamps
        logger.info("Sorting sweeps based on timestamp")
        sweeps = sorted(sweeps, key=lambda x: x.get_scan().get_epochs()[0])

        # Count the number of datasets from each day
        from collections import Counter

        counter = Counter()
        for s in sweeps:
            timestamp = s.get_scan().get_epochs()[0]
            timestamp = datetime.datetime.fromtimestamp(timestamp)
            timestamp = timestamp.strftime("%Y-%m-%d")
            counter[timestamp] += 1

        # Print the number of datasets on each day
        for timestamp in sorted(counter.keys()):
            logger.info("%d datasets collected on %s" % (counter[timestamp], timestamp))

        # Loop though and see if any models might be shared
        b_list = [s.get_beam() for s in sweeps]
        d_list = [s.get_detector() for s in sweeps]
        g_list = [s.get_goniometer() for s in sweeps]
        b_index = []
        d_index = []
        g_index = []
        for i in range(len(sweeps)):
            b = b_list[i]
            d = d_list[i]
            g = g_list[i]
            bn = i
            dn = i
            gn = i
            if i > 0:
                bj = b_index[-1]
                dj = d_index[-1]
                gj = g_index[-1]
                if b.is_similar_to(b_list[bj]):
                    bn = bj
                if d.is_similar_to(d_list[dj]):
                    dn = dj
                if g.is_similar_to(g_list[gj]):
                    gn = gj
            b_index.append(bn)
            d_index.append(dn)
            g_index.append(gn)

        # Print a table of possibly shared models
        from libtbx.table_utils import format as table

        rows = [["Sweep", "ID", "Beam", "Detector", "Goniometer", "Date", "Time"]]
        for i in range(len(sweeps)):
            timestamp = sweeps[i].get_scan().get_epochs()[0]
            timestamp = datetime.datetime.fromtimestamp(timestamp)
            date_str = timestamp.strftime("%Y-%m-%d")
            time_str = timestamp.strftime("%H:%M:%S")
            row = [
                "%s" % sweeps[i].get_template(),
                "%s" % i,
                "%s" % b_index[i],
                "%s" % d_index[i],
                "%s" % g_index[i],
                "%s" % date_str,
                "%s" % time_str,
            ]
            rows.append(row)
        logger.info(table(rows, has_header=True, justify="left", prefix=" "))

        # Print the time
        logger.info("Time Taken: %f" % (time() - start_time))


if __name__ == "__main__":
    from dials.util import halraiser

    try:
        script = Script()
        script.run()
    except Exception as e:
        halraiser(e)

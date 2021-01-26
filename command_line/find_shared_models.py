import datetime
import logging
from collections import Counter

from dxtbx.imageset import ImageSequence
from libtbx.phil import parse

from dials.util import log, show_mail_handle_errors, tabulate
from dials.util.options import OptionParser, flatten_experiments
from dials.util.version import dials_version

logger = logging.getLogger("dials.command_line.find_shared_models")

help_message = """

This program attempts to find sets of images with shared models

Examples::

  dials.find_shared_models models.expt

  dials.find_shared_models experiments1.expt experiments2.expt
"""

# Set the phil scope
phil_scope = parse(
    """
  output {
    log = 'dials.find_shared_models.log'
      .type = str
      .help = "The log filename"
  }
""",
    process_includes=True,
)


class Script:
    """A class for running the script."""

    def __init__(self):
        """Initialise the script."""
        # The script usage
        usage = "dials.find_shared_models [options] [param.phil] models.expt"

        # Initialise the base class
        self.parser = OptionParser(
            usage=usage,
            phil=phil_scope,
            epilog=help_message,
            read_experiments=True,
            check_format=False,
        )

    def run(self, args=None):
        """Execute the script."""

        # Parse the command line
        params, options = self.parser.parse_args(args, show_diff_phil=False)

        # Configure the logging
        log.config(verbosity=options.verbose, logfile=params.output.log)
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

        # Get the list of sequences
        sequences = []
        for experiment in experiments:
            if isinstance(experiment.imageset, ImageSequence):
                sequences.append(experiment.imageset)
        logger.info("Number of sequences = %d", len(sequences))

        # Sort the sequences by timestamps
        logger.info("Sorting sequences based on timestamp")
        sequences = sorted(sequences, key=lambda x: x.get_scan().get_epochs()[0])

        # Count the number of datasets from each day
        counter = Counter()
        for s in sequences:
            timestamp = s.get_scan().get_epochs()[0]
            timestamp = datetime.datetime.fromtimestamp(timestamp)
            timestamp = timestamp.strftime("%Y-%m-%d")
            counter[timestamp] += 1

        # Print the number of datasets on each day
        for timestamp in sorted(counter.keys()):
            logger.info("%d datasets collected on %s", counter[timestamp], timestamp)

        # Loop though and see if any models might be shared
        b_list = [s.get_beam() for s in sequences]
        d_list = [s.get_detector() for s in sequences]
        g_list = [s.get_goniometer() for s in sequences]
        b_index = []
        d_index = []
        g_index = []
        for i in range(len(sequences)):
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
        rows = [["Sequence", "ID", "Beam", "Detector", "Goniometer", "Date", "Time"]]
        for i in range(len(sequences)):
            timestamp = sequences[i].get_scan().get_epochs()[0]
            timestamp = datetime.datetime.fromtimestamp(timestamp)
            date_str = timestamp.strftime("%Y-%m-%d")
            time_str = timestamp.strftime("%H:%M:%S")
            row = [
                "%s" % sequences[i].get_template(),
                "%s" % i,
                "%s" % b_index[i],
                "%s" % d_index[i],
                "%s" % g_index[i],
                "%s" % date_str,
                "%s" % time_str,
            ]
            rows.append(row)
        logger.info(tabulate(rows, headers="firstrow"))


@show_mail_handle_errors()
def run(args=None):
    script = Script()
    script.run(args)


if __name__ == "__main__":
    run()

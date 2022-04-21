from __future__ import annotations

import logging
import sys

from libtbx.phil import parse

from dials.algorithms.statistics.cc_half_algorithm import CCHalfFromDials, CCHalfFromMTZ
from dials.util import log, show_mail_handle_errors
from dials.util.options import ArgumentParser, reflections_and_experiments_from_files

logger = logging.getLogger("dials.command_line.compute_delta_cchalf")

help_message = """

This program computes the delta cchalf excluding images
"""

# Set the phil scope
phil_scope = parse(
    """

  input {

    mtzfile = None
      .type = str
      .help = "We can also import an MTZ file"

  }

  mode = *dataset image_group
    .type = choice
    .help = "Perform analysis on whole datasets or batch groups"

  group_size = 10
    .type = int(value_min=1)
    .help = "The number of images to group together when calculating delta"
            "cchalf in image_group mode"

  mtz {
    batch_offset = None
      .type = int
      .help = "The batch offset between consecutive datasets in the mtz file."
  }

  output {

    experiments = "filtered.expt"
      .type = str
      .help = "The filtered experiments file"

    reflections = "filtered.refl"
      .type = str
      .help = "The filtered reflections file"

    table = "delta_cchalf.dat"
      .type = str
      .help = "A file with delta cchalf values"
    html = "compute_delta_cchalf.html"
      .type = str
      .help = "HTML filename for report of results."
  }

  nbins = 10
    .type = int(value_min=1)
    .help = "The number of resolution bins to use"

  dmin = None
    .type = float
    .help = "The maximum resolution"

  dmax = None
    .type = float
    .help = "The minimum resolution"

  stdcutoff = 4.0
    .type = float
    .help = "Datasets with a ΔCC½ below (mean - stdcutoff*std) are removed"

  output {
    log = 'dials.compute_delta_cchalf.log'
      .type = str
      .help = "The log filename"
  }
"""
)


@show_mail_handle_errors()
def run(args=None, phil=phil_scope):
    """Run the command-line script."""

    usage = "dials.compute_delta_cchalf [options] scaled.expt scaled.refl"

    parser = ArgumentParser(
        usage=usage,
        phil=phil,
        epilog=help_message,
        read_experiments=True,
        read_reflections=True,
        check_format=False,
    )

    params, _ = parser.parse_args(args=args, show_diff_phil=False)

    log.config(logfile=params.output.log)

    reflections, experiments = reflections_and_experiments_from_files(
        params.input.reflections, params.input.experiments
    )

    if not experiments and not reflections:
        if not params.input.mtzfile:
            parser.print_help()
            return
        else:
            try:
                script = CCHalfFromMTZ(params, params.input.mtzfile)
            except ValueError as e:
                sys.exit(f"Error: {e}")
    else:
        if not experiments or not reflections:
            parser.print_help()
            return
        else:
            if not len(reflections) == 1:
                exit("Only one reflection table can be provided")
            n_datasets = len(set(reflections[0]["id"]).difference({-1}))
            if n_datasets != len(experiments):
                exit(
                    """
The number of experiments (%s) does not match the number
of datasets in the reflection table (%s)
""",
                    len(experiments),
                    n_datasets,
                )
            try:
                script = CCHalfFromDials(params, experiments, reflections[0])
            except ValueError as e:
                sys.exit(f"Error: {e}")

    script.run()
    script.output()


if __name__ == "__main__":
    run()

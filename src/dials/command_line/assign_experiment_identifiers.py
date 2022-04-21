from __future__ import annotations

import sys

from libtbx import phil

from dials.array_family import flex
from dials.util import Sorry, show_mail_handle_errors
from dials.util.multi_dataset_handling import (
    assign_unique_identifiers,
    parse_multiple_datasets,
)
from dials.util.options import ArgumentParser, reflections_and_experiments_from_files

help_message = """Command line script which assigns experiment identifiers
to reflections and experiments and saves them back to disk.
"""

phil_scope = phil.parse(
    """
  identifiers = None
    .type = strings
    .help = "User specified identifiers to use, must be a list of strings equal"
            "to the number of datasets."
  output {
    reflections = assigned.refl
      .type = str
    experiments = assigned.expt
      .type = str
  }
"""
)


@show_mail_handle_errors()
def run(args=None):
    """Run assign experiment identifiers from the command line."""
    usage = (
        """Usage: dials.assign_experiment_identifiers observations.refl models.expt"""
    )
    parser = ArgumentParser(
        usage=usage,
        read_experiments=True,
        read_reflections=True,
        phil=phil_scope,
        check_format=False,
        epilog=help_message,
    )
    params, _ = parser.parse_args(args=args, show_diff_phil=False)

    if not params.input.experiments or not params.input.reflections:
        parser.print_help()
        sys.exit()

    reflections, experiments = reflections_and_experiments_from_files(
        params.input.reflections, params.input.experiments
    )

    reflections = parse_multiple_datasets(reflections)
    if len(experiments) != len(reflections):
        raise Sorry(
            "Mismatched number of experiments and reflection tables found: %s & %s."
            % (len(experiments), len(reflections))
        )
    try:
        experiments, reflections = assign_unique_identifiers(
            experiments, reflections, params.identifiers
        )
    except ValueError as e:
        raise Sorry(e)
    print(f"assigned identifiers: {list(experiments.identifiers())}")

    experiments.as_file(params.output.experiments)
    joint_table = flex.reflection_table()
    for reflection_table in reflections:
        joint_table.extend(reflection_table)
    joint_table.as_file(params.output.reflections)


if __name__ == "__main__":
    run()

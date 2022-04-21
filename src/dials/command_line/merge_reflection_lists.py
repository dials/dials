# LIBTBX_SET_DISPATCHER_NAME dials.merge_reflection_lists


from __future__ import annotations

import sys

from libtbx.phil import parse

from dials.util import show_mail_handle_errors
from dials.util.command_line import Command
from dials.util.options import ArgumentParser

# Create the help message
help_message = """

  This program allows you to combine reflection lists. Note that merge in this
  case refers to merging lists rather than merging equivalent reflections.

  Reflection tables can either be extended (i.e. given two reflection tables,
  "A" and "B", the result would have length = len("A") + len("B")). Or values
  in one table can be overwritten with values from the other table. Any columns
  present in "A" but not present in "B" (or vice versa) will be preserved.
"""

# Create the phil parameters
phil_scope = parse(
    """

  output = merged.refl
    .type = str
    .help = "The output file"

  method = *update extend
    .type = choice
    .help = "The method of merging"
"""
)


class Script:
    """A class to encapsulate the script."""

    def __init__(self):
        """Initialise the script."""
        # The script usage
        usage = "usage: dials.merge_reflection_lists [options] /path/to/image/reflection/files"
        self.parser = ArgumentParser(
            epilog=help_message, usage=usage, phil=phil_scope, read_reflections=True
        )

    def run(self, args=None):
        """Run the script."""
        # Parse the command line arguments
        params, options = self.parser.parse_args(args, show_diff_phil=True)
        if len(params.input.reflections) == 0:
            self.parser.print_help()
            return
        if len(params.input.reflections) <= 1:
            sys.exit("more than 1 reflection table must be specified")
        tables = [p.data for p in params.input.reflections]

        # Get the number of rows and columns
        nrows = [t.nrows() for t in tables]
        ncols = [t.ncols() for t in tables]

        # Merge the reflection lists
        if params.method == "update":
            assert all(n == nrows[0] for n in nrows[1:])
            table = tables[0]
            for t in tables[1:]:
                table.update(t)
        elif params.method == "extend":
            assert all(n == ncols[0] for n in ncols[1:])
            table = tables[0]
            for t in tables[1:]:
                table.extend(t)
        else:
            raise RuntimeError(f"unknown method, {params.method}")

        # Write the reflections to the file
        Command.start(f"Writing {len(table)} reflections to {params.output}")
        table.as_file(params.output)
        Command.end(f"Wrote {len(table)} reflections to {params.output}")


@show_mail_handle_errors()
def run(args=None):
    script = Script()
    script.run(args)


if __name__ == "__main__":
    run()

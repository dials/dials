#!/usr/bin/env python
#
# merge_reflection_lists.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

# LIBTBX_SET_DISPATCHER_NAME dev.dials.merge_reflection_lists

from __future__ import absolute_import, division, print_function

from libtbx.phil import parse

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


class Script(object):
    """ A class to encapsulate the script. """

    def __init__(self):
        """ Initialise the script. """
        from dials.util.options import OptionParser
        import libtbx.load_env

        # The script usage
        usage = (
            "usage: %s [options] /path/to/image/reflection/files"
            % libtbx.env.dispatcher_name
        )
        self.parser = OptionParser(
            epilog=help_message, usage=usage, phil=phil_scope, read_reflections=True
        )

    def run(self):
        """ Run the script. """
        from dials.util.command_line import Command
        from dials.util import Sorry

        # Parse the command line arguments
        params, options = self.parser.parse_args(show_diff_phil=True)
        if len(params.input.reflections) == 0:
            self.parser.print_help()
            return
        if len(params.input.reflections) <= 1:
            raise Sorry("more than 1 reflection table must be specified")
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
            raise RuntimeError("unknown method, %s" % params.method)

        # Write the reflections to the file
        Command.start("Writing %d reflections to %s" % (len(table), params.output))
        table.as_pickle(params.output)
        Command.end("Wrote %d reflections to %s" % (len(table), params.output))


if __name__ == "__main__":
    from dials.util import halraiser

    try:
        script = Script()
        script.run()
    except Exception as e:
        halraiser(e)

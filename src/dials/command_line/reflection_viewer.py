# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1

from __future__ import annotations

import dials.util

help_message = """

This program is used to view the reflections with debugging purposes.
This program does not perform any calculation ... just visualizations

Example for invoking from CLI:

dials.reflection_viewer observations.refl
"""


class Script:
    """The debugging visualization program."""

    def __init__(self):
        """Initialise the script."""
        from dials.util.options import ArgumentParser

        usage = "dials.reflection_viewer [options] experiment.expt"

        # Create the parser
        self.parser = ArgumentParser(
            usage=usage, epilog=help_message, read_reflections=True
        )

    def run(self, args=None):

        from dials.util.options import flatten_reflections
        from dials.viewer.viewer_interface import extract_n_show

        # Parse the command line
        params, options = self.parser.parse_args(args, show_diff_phil=True)
        table = flatten_reflections(params.input.reflections)
        if len(table) == 0:
            self.parser.print_help()
            return

        extract_n_show(table[0])


@dials.util.show_mail_handle_errors()
def run(args=None):
    script = Script()
    script.run(args)


if __name__ == "__main__":
    run()

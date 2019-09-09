from __future__ import absolute_import, division, print_function

# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1


help_message = """

This program is used to view the reflections with debugging purposes.
This program does not perform any calculation ... just visualizations

Example for invoking from CLI:

dials.reflection_viewer observations.refl
"""


class Script(object):
    """ The debugging visualization program. """

    def __init__(self):
        """Initialise the script."""
        from dials.util.options import OptionParser

        usage = "dials.reflection_viewer [options] experiment.expt"

        # Create the parser
        self.parser = OptionParser(
            usage=usage, epilog=help_message, read_reflections=True
        )

    def run(self):

        from dials.util.options import flatten_reflections
        from dials.viewer.viewer_interface import extract_n_show

        # Parse the command line
        params, options = self.parser.parse_args(show_diff_phil=True)
        table = flatten_reflections(params.input.reflections)
        if len(table) == 0:
            self.parser.print_help()
            return

        extract_n_show(table[0])


if __name__ == "__main__":
    script = Script()
    script.run()

#!/usr/bin/env python
#
# dials.reflection_viewer.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst and Luis Fuentes-Montero (Luiso)
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import absolute_import, division, print_function

# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1


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
        import libtbx.load_env

        # The script usage
        usage = "usage: %s [options] experiment.expt" % libtbx.env.dispatcher_name

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

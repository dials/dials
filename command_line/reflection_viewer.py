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

from __future__ import division

help_message = '''

This program is used to view the reflections with debugging purposes.
This program does not perform any calculation ... just visualizations

Example for invoking from CLI:

dials.reflection_viewer My_Reflections.pickle

'''

class Script(object):
  ''' The debugging visualization program. '''

  def __init__(self):
    '''Initialise the script.'''

    # The script usage
    usage  = "usage: %prog [options] experiment.json"

    # Create the parser
    self.parser = OptionParser(
      usage=usage,
      epilog=help_message)

  def run(self):

    ''' Show the reflections one by one in an interactive way '''
    from dials.util.command_line import Importer
    from dials.viewer.reflection_view import viewer_App
    from dials.array_family import flex

    # Parse the command line
    params, options, args = self.parser.parse_args(show_diff_phil=True)
    if len(args) == 0:
      self.parser.print_help()
      return

    importer = Importer(args, include=["reflections", "experiments"])
    my_tables = importer.reflections

    ''' opens and closes the viewer for each new reflection table '''
    for table in my_tables:
      My_app = viewer_App(redirect=False)
      My_app.table_in(table)
      My_app.MainLoop()
      My_app.Destroy()

if __name__ == '__main__':
  script = Script()
  script.run()

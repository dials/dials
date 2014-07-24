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
from dials.util.script import ScriptRunner

help_message = '''

This program is used to view the reflections with debugging purposes.
This program does not perform any calculation ... just visualizations

Example for invoking from CLI:

dials.reflection_viewer My_Reflections.pickle

'''

class Script(ScriptRunner):
  ''' The debugging visualization program. '''

  def __init__(self):
    '''Initialise the script.'''

    # The script usage
    usage  = "usage: %prog [options] experiment.json"

    # Initialise the base class
    ScriptRunner.__init__(self, usage=usage,
                          epilog=help_message,
                          home_scope="integration")

  def main(self, params, options, args):

    ''' Show the reflections one by one in an interactive way '''
    from dials.util.command_line import Importer
    from dials.viewer.reflection_view import viewer_App
    from dials.array_family import flex
    if len(args) == 0:
      self.config().print_help()
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

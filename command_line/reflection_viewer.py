#!/usr/bin/env python
#
# dials.integrate.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division
from dials.util.script import ScriptRunner

help_message = '''

This program is used to view the reflections on the diffraction images. It
is called with an experiment list
Examples:

TO DO

'''

class Script(ScriptRunner):
  ''' The integration program. '''

  def __init__(self):
    '''Initialise the script.'''

    # The script usage
    usage  = "usage: %prog [options] experiment.json"

    # Initialise the base class
    ScriptRunner.__init__(self, usage=usage,
                          epilog=help_message,
                          home_scope="integration")

  def main(self, params, options, args):
    ''' Perform the integration. '''
    from dials.util.command_line import Importer
    from dials.viewer.reflection_view import viewer_App
    from dials.array_family import flex


    importer = Importer(args, include=["reflections", "experiments"])
    print importer.unhandled_arguments
    print importer.experiments
    print importer.reflections

    my_tables = importer.reflections
    print "len =", len(my_tables)

    for table in my_tables:
      My_app = viewer_App(redirect=False)
      My_app.table_in(table)
      My_app.MainLoop()
      My_app.Destroy()

if __name__ == '__main__':
  script = Script()
  script.run()

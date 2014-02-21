#!/usr/bin/env python
#
# dials.process.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division
from dials.util.script import ScriptRunner


class Script(ScriptRunner):
  '''A class for running the script.'''

  def __init__(self):
    '''Initialise the script.'''

    # The script usage
    usage = "usage: %prog [options] [param.phil] datablock.json"

    # Initialise the base class
    ScriptRunner.__init__(self, usage=usage)

    # Add a verbosity option
    self.config().add_option(
        "-v",
        action="count", default=1,
        help="set verbosity level; -vv gives verbosity level 2")

  def main(self, params, options, args):
    '''Execute the script.'''
    from dials.util.command_line import Importer
    from dials.util.command_line import Command
    from time import time
    st = time()

    # Preamble stuff
    print '*' * 80
    print ''
    print '                       mmmm   mmmmm    mm   m       mmmm            '
    print '                       #   "m   #      ##   #      #"   "           '
    print '                      m#mm  #   #     #  #  #      "#mmm            '
    print '                       #    #   #     #mm#  #          "#           '
    print '                       #mmm"  mm#mm  #    # #mmmmm "mmm#"           '
    print ''
    print 'Launching dials.process'
    print ''
    print 'The following tasks will be performed:'
    print ' 1) Strong spots will be found (dials.find_spots)'
    print ' 2) The strong spots will be indexed (dials.index)'
    print ' 3) A profile model will be created (dials.create_profile_model)'
    print ' 4) The reflections will be integrated (dials.integrate)'
    print ''
    print 'Please be patient, this may take a few minutes'
    print ''
    print '*' * 80
    print ''

    # Import stuff
    Command.start('Importing datablocks')
    importer = Importer(args)
    assert(len(importer.datablocks) == 1)
    datablock = importer.datablocks[0]
    Command.end('Imported datablocks')
    print ''

    # Find the strong spots
    observed = self.find_spots(datablock)

    # Index the strong spots
    experiments, indexed = self.index(datablock, observed)

    # Create the profile model
    profile = self.create_profile_model(experiments, indexed)

    # Integrate the reflections
    integrated = self.integrate(experiments, profile, indexed)

    # Total Time
    print ""
    print "Total Time Taken = %f seconds" % (time() - st)

  def find_spots(self, datablock):
    from time import time
    st = time()

    print '*' * 80
    print 'Finding Strong Spots'
    print '*' * 80

    print ''
    print 'Time Taken = %f seconds' % (time() - st)
    return None

  def index(self, datablock, observed):
    from time import time
    st = time()

    print '*' * 80
    print 'Indexing Strong Spots'
    print '*' * 80

    print ''
    print 'Time Taken = %f seconds' % (time() - st)
    return None, None

  def create_profile_model(self, experiments, indexed):
    from time import time
    st = time()

    print '*' * 80
    print 'Creating Profile Model'
    print '*' * 80

    print ''
    print 'Time Taken = %f seconds' % (time() - st)
    return None

  def integrate(self, experiments, profile, indexed):
    from time import time
    st = time()

    print '*' * 80
    print 'Integrating Reflections'
    print '*' * 80

    print ''
    print 'Time Taken = %f seconds' % (time() - st)
    return None


if __name__ == '__main__':
  script = Script()
  script.run()

#!/usr/bin/env python
#
# plot_scan_varying_crystal.py
#
#  Copyright (C) 2014 Diamond Light Source
#
#  Author: David G. Waterman
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division


class ScriptRunner(object):
  '''Class to run script.'''

  def __init__(self, crystals):
    '''Setup the script.'''

    # Filename data
    self.crystals = crystals

  def __call__(self):
    '''Run the script.'''
    #from dials.model.data import ReflectionList # import dependency
    #from dials.util.command_line import Command

    for icrystal, crystal in enumerate(self.crystals):

      if crystal.num_scan_points == 0:
        print "Ignoring scan-static crystal"
        continue

      print "Image\ta\tb\tc\talpha\tbeta\tgamma"
      msg = "\t".join(["%.3f"] * 7)
      for t in xrange(crystal.num_scan_points):
        uc = crystal.get_unit_cell_at_scan_point(t)
        print msg % ((t,) + uc.parameters())

    print "TODO: misset angles around user-supplied axes"

if __name__ == '__main__':
  import sys
  from dials.util.command_line import Importer
  args = sys.argv[1:]
  importer = Importer(args, check_format=False)
  try:
    crystals = importer.experiments.crystals()
  except AttributeError:
    print "No crystals found in the input"
    raise

  assert len(importer.unhandled_arguments) == 0

  runner = ScriptRunner(
      crystals=crystals)

  # Run the script
  runner()

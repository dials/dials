#!/usr/bin/env python
#
# dials.simulate.py
#
#  Copyright (C) 2013 Diamond Light Source, Graeme Winter
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division
from dials.util.script import ScriptRunner

class Script(ScriptRunner):

  def __init__(self):

    usage  = "usage: %prog [options] [param.phil] " \
             "sweep.json crystal.json intensities.mtz"

    ScriptRunner.__init__(self, usage=usage)

    self.config().add_option(
        '--output',
        dest = 'output',
        type = 'string', default = 'simulation.pickle',
        help = 'Set the filename for simulated reflection file.')

    self.config().add_option(
        "-v", "--verbosity",
        action="count", default=0,
        help="set verbosity level; -vv gives verbosity level 2.")

  def main(self, params, options, args):
    # FIXME import simulation code
    from dials.model.serialize import load, dump
    from iotbx import mtz
    import cPickle as pickle
    from dials.util.command_line import Importer

    importer = Importer(args)
    if len(importer.imagesets) == 0 and len(importer.crystals) == 0:
      self.config().print_help()
      return
    if len(importer.imagesets) != 1:
      raise RuntimeError('need 1 sweep: %d given' % len(importer.imagesets))
    if len(importer.crystals) != 1:
      raise RuntimeError('need 1 crystal: %d given' % len(importer.crystals))
    sweep = importer.imagesets[0]
    crystal = importer.crystals[0]
    data = mtz.object(args[2])

    goniometer = sweep.get_goniometer()
    detector = sweep.get_detector()
    beam = sweep.get_beam()
    scan = sweep.get_scan()

    # FIXME generate predictions for requested reflections => generate a
    # reflection list

    # FIXME calculate shoebox sizes: take parameters from params & transform
    # from reciprocal space to image space to decide how big a shoe box to use

    # FIXME now for each reflection perform the simulation

    # FIXME now for each reflection add background


    return

if __name__ == '__main__':
  script = Script()
  script.run()

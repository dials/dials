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
    from dials.algorithms.integration import ReflectionPredictor

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

    predict = ReflectionPredictor()
    predicted = predict(sweep, crystal)

    from dials.scratch.jmp.container.reflection_table import ReflectionTable
    
    # FIXME calculate shoebox sizes: take parameters from params & transform
    # from reciprocal space to image space to decide how big a shoe box to use

    table = ReflectionTable()
    table['miller_index'] = predicted.miller_index()
    indexer = table.index_map('miller_index')

    candidates = []

    unique = sorted(indexer)
    
    for h, k, l in unique:

      try:
        for _h in h - 1, h + 1:
          if not indexer[(_h, k, l)]:
            raise ValueError, 'missing'
        for _k in k - 1, k + 1:
          if not indexer[(h, _k, l)]:
            raise ValueError, 'missing'
        for _l in l - 1, l + 1:
          if not indexer[(h, k, _l)]:
            raise ValueError, 'missing'
        candidates.append((h, k, l))
      except ValueError, e:
        continue

    from dials.algorithms.simulation.utils import build_prediction_matrix

    # FIXME should this not be handled by magic by the ScriptRunner?
    
    from dials.algorithms.simulation.generate_test_reflections import \
     master_phil
    from libtbx.phil import command_line
    cmd = command_line.argument_interpreter(master_params=master_phil)
    working_phil = cmd.process_and_fetch(args=args[3:])
    params = working_phil.extract()
    
    node_size = params.rs_node_size
    window_size = params.rs_window_size
    
    from dials.model.data import ReflectionList

    useful = ReflectionList()
    d_matrices = []

    for h, k, l in candidates:      
      hkl = predicted[indexer[(h, k, l)][0]]
      _x = hkl.image_coord_px[0]
      _y = hkl.image_coord_px[1]
      _z = hkl.frame_number

      # build prediction matrix
      mhkl = predicted[indexer[(h - 1, k, l)][0]]
      phkl = predicted[indexer[(h + 1, k, l)][0]]
      hmkl = predicted[indexer[(h, k - 1, l)][0]]
      hpkl = predicted[indexer[(h, k + 1, l)][0]]
      hkml = predicted[indexer[(h, k, l - 1)][0]]
      hkpl = predicted[indexer[(h, k, l + 1)][0]]
      d = build_prediction_matrix(hkl, mhkl, phkl, hmkl, hpkl, hkml, hkpl)
      d_matrices.append(d)
      
      # construct the shoebox parameters
      x, y, z = [], [], []

      for dh in (1, 0, 0), (0, 1, 0), (0, 0, 1):
        dxyz = -1 * window_size * d * dh
        x.append(dxyz[0] + _x)
        y.append(dxyz[1] + _y)
        z.append(dxyz[2] + _z)
        dxyz = window_size * d * dh
        x.append(dxyz[0] + _x)
        y.append(dxyz[1] + _y)
        z.append(dxyz[2] + _z)

      import math
        
      hkl.bounding_box = (int(math.floor(min(x))), int(math.floor(max(x)) + 1),
                          int(math.floor(min(y))), int(math.floor(max(y)) + 1),
                          int(math.floor(min(z))), int(math.floor(max(z)) + 1))
      useful.append(hkl)

    from dials.algorithms import shoebox
    shoebox.allocate(useful)

    useful = useful[:10]
    
    import random
    
    # FIXME now for each reflection perform the simulation
    for j, refl in enumerate(useful):
      d = d_matrices[j]
      for j in range(params.counts):
        dh = random.gauss(0, node_size)
        dk = random.gauss(0, node_size)
        dl = random.gauss(0, node_size)
        _x = refl.image_coord_px[0]
        _y = refl.image_coord_px[1]
        _z = refl.frame_number
        xyz = d * (dh, dk, dl)
        x = int(math.floor(xyz[0] + _x))
        y = int(math.floor(xyz[1] + _y))
        z = int(math.floor(xyz[2] + _z))
        if x < refl.bounding_box[0] or x >= refl.bounding_box[1]:
          continue
        if y < refl.bounding_box[2] or y >= refl.bounding_box[3]:
          continue
        if z < refl.bounding_box[4] or z >= refl.bounding_box[5]:
          continue
        x -= refl.bounding_box[0]
        y -= refl.bounding_box[2]
        z -= refl.bounding_box[4]
        refl.shoebox[(z, y, x)] += 1

    pickle.dump(useful, open('useful.pickle', 'w'))
        
    # FIXME now for each reflection add background


    return

if __name__ == '__main__':
  script = Script()
  script.run()

#!/usr/bin/env python
#
# dials.simulate.py
#
#  Copyright (C) 2013 Diamond Light Source, Graeme Winter
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import absolute_import, division, print_function

# LIBTBX_SET_DISPATCHER_NAME dev.dials.simulate

class Script(object):

  def __init__(self):
    from dials.util.options import OptionParser
    from libtbx.phil import parse

    usage  = "usage: %prog [options] [param.phil] " \
             "sweep.json crystal.json intensities.mtz"

    phil_scope = parse('''
      output = simulated.pickle
        .type = str
        .help = "The filename for the simulated reflections"
    ''')

    # Create the parser
    self.parser = OptionParser(
      usage=usage,
      phil=self.phil_scope())

  @staticmethod
  def map_to_image_space(refl, d, dhs, dks, dls):
    from scitbx.array_family import flex
    d_elems = d.elems
    bb = refl.bounding_box
    dxs = d_elems[0] * dhs + d_elems[1] * dks + d_elems[2] * dls
    dys = d_elems[3] * dhs + d_elems[4] * dks + d_elems[5] * dls
    dzs = d_elems[6] * dhs + d_elems[7] * dks + d_elems[8] * dls
    xs = flex.floor(dxs + refl.image_coord_px[0]).iround() - bb[0]
    ys = flex.floor(dys + refl.image_coord_px[1]).iround() - bb[2]
    zs = flex.floor(dzs + refl.frame_number).iround() - bb[4]
    xyz = flex.vec3_int(zs, ys, xs)
    xyz = xyz.select((xs >= 0 and xs < (bb[1] - bb[0])) &
                     (ys >= 0 and ys < (bb[3] - bb[2])) &
                     (zs >= 0 and zs < (bb[5] - bb[4])))
    for _xyz in xyz:
      refl.shoebox[_xyz] += 1

    return

  def main(self):
    # FIXME import simulation code
    import six.moves.cPickle as pickle
    import math
    from dials.util.command_line import Importer
    from dials.algorithms.integration import ReflectionPredictor
    from libtbx.utils import Sorry

    # Parse the command line
    params, options, args = self.parser.parse_args()

    importer = Importer(args)
    if len(importer.imagesets) == 0 and len(importer.crystals) == 0:
      self.config().print_help()
      return
    if len(importer.imagesets) != 1:
      raise Sorry('need 1 sweep: %d given' % len(importer.imagesets))
    if len(importer.crystals) != 1:
      raise Sorry('need 1 crystal: %d given' % len(importer.crystals))
    sweep = importer.imagesets[0]
    crystal = importer.crystals[0]

    # generate predictions for possible reflections => generate a
    # reflection list

    predict = ReflectionPredictor()
    predicted = predict(sweep, crystal)

    # sort with James's reflection table: should this not go somewhere central?
    from dials.scratch.jmp.container.reflection_table import ReflectionTable

    # calculate shoebox sizes: take parameters from params & transform
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
            raise ValueError('missing')
        for _k in k - 1, k + 1:
          if not indexer[(h, _k, l)]:
            raise ValueError('missing')
        for _l in l - 1, l + 1:
          if not indexer[(h, k, _l)]:
            raise ValueError('missing')
        candidates.append((h, k, l))
      except ValueError:
        continue

    from dials.algorithms.simulation.utils import build_prediction_matrix

    from dials.algorithms.simulation.generate_test_reflections import \
     master_phil
    from libtbx.phil import command_line
    cmd = command_line.argument_interpreter(master_params=master_phil)
    working_phil = cmd.process_and_fetch(args=args[2:])
    params = working_phil.extract()

    node_size = params.rs_node_size
    window_size = params.rs_window_size
    reference = params.integrated_data_file
    scale = params.integrated_data_file_scale

    if reference:
      counts_database = { }
      from iotbx import mtz
      m = mtz.object(reference)
      mi = m.extract_miller_indices()
      i = m.extract_reals('IMEAN').data
      s = m.space_group().build_derived_point_group()
      for j in range(len(mi)):
        for op in s.all_ops():
          hkl = tuple(map(int, op * mi[j]))
          counts = max(0, int(math.floor(i[j] * scale)))
          counts_database[hkl] = counts
          counts_database[(-hkl[0], -hkl[1], -hkl[2])] = counts
    else:
      def constant_factory(value):
        import itertools
        return itertools.repeat(value).next
      from collections import defaultdict
      counts_database = defaultdict(constant_factory(params.counts))

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

      # construct the shoebox parameters: outline the ellipsoid
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

      hkl.bounding_box = (int(math.floor(min(x))), int(math.floor(max(x)) + 1),
                          int(math.floor(min(y))), int(math.floor(max(y)) + 1),
                          int(math.floor(min(z))), int(math.floor(max(z)) + 1))
      try:
        counts = counts_database[hkl.miller_index]
        useful.append(hkl)
      except KeyError:
        continue

    from dials.algorithms import shoebox
    shoebox.allocate(useful)

    from dials.util.command_line import ProgressBar
    p = ProgressBar(title = 'Generating shoeboxes')

    # now for each reflection perform the simulation
    for j, refl in enumerate(useful):
      p.update(j * 100.0 / len(useful))
      d = d_matrices[j]

      from scitbx.random import variate, normal_distribution
      g = variate(normal_distribution(mean = 0, sigma = node_size))
      counts = counts_database[refl.miller_index]
      dhs = g(counts)
      dks = g(counts)
      dls = g(counts)
      self.map_to_image_space(refl, d, dhs, dks, dls)

    p.finished('Generated %d shoeboxes' % len(useful))

    # now for each reflection add background
    from dials.algorithms.simulation.generate_test_reflections import \
     random_background_plane

    p = ProgressBar(title = 'Generating background')
    for j, refl in enumerate(useful):
      p.update(j * 100.0 / len(useful))
      if params.background:
        random_background_plane(refl.shoebox, params.background, 0.0, 0.0, 0.0)
      else:
        random_background_plane(
          refl.shoebox, params.background_a, params.background_b,
          params.background_c, params.background_d)

    p.finished('Generated %d backgrounds' % len(useful))
    if params.output.all:
      with open(params.output.all, 'wb') as fh:
        pickle.dump(useful, fh, pickle.HIGHEST_PROTOCOL)

if __name__ == '__main__':
  from dials.util import halraiser
  try:
    script = Script()
    script.run()
  except Exception as e:
    halraiser(e)

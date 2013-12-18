#!/usr/bin/env python
#
# export_xds.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: Richard Gildea
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import division


def run(args):
  from libtbx.phil import command_line
  from dials.util.command_line import Importer

  importer = Importer(args)
  imagesets = importer.imagesets
  crystals = importer.crystals
  assert len(imagesets) == len(crystals)
  reflection_lists = importer.reflections
  args = importer.unhandled_arguments

  from dxtbx.serialize import xds
  from iotbx.xds import spot_xds

  for i in range(len(imagesets)):
    suffix = ""
    if len(imagesets) > 1:
      suffix = "_%i" %(i+1)

    imageset = imagesets[i]
    crystal_model = crystals[i]
    crystal_model = crystal_model.change_basis(
      crystal_model.get_space_group().info()\
        .change_of_basis_op_to_reference_setting())
    A = crystal_model.get_A()
    A_inv = A.inverse()
    real_space_a = A_inv.elems[:3]
    real_space_b = A_inv.elems[3:6]
    real_space_c = A_inv.elems[6:9]
    to_xds = xds.to_xds(imageset)
    with open('XDS%s.INP' %suffix, 'wb') as f:
      to_xds.XDS_INP(out=f, job_card="XYCORR INIT DEFPIX INTEGRATE CORRECT")
    with open('XPARM%s.XDS' %suffix, 'wb') as f:
      to_xds.xparm_xds(
        real_space_a, real_space_b, real_space_c,
        crystal_model.get_space_group().type().number(),
        out=f)

  for i, reflections in enumerate(reflection_lists):
    suffix = ""
    if len(reflection_lists) > 1:
      suffix = "_%i" %(i+1)

    centroids = []
    intensities = []
    miller_indices = []
    miller_indices_excluding_zero = []

    for refl in reflections:
      centroids.append(refl.centroid_position)
      intensities.append(refl.intensity)
      miller_indices_excluding_zero.append(refl.miller_index)
      if refl.miller_index != (0,0,0):
        miller_indices.append(refl.miller_index)

    if len(miller_indices) == 0:
      miller_indices = None
    xds_writer = spot_xds.writer(centroids=centroids,
                                 intensities=intensities,
                                 miller_indices=miller_indices)
    xds_writer.write_file(filename='SPOTS%s.XDS' %suffix)


if __name__ == '__main__':
  import sys
  run(sys.argv[1:])

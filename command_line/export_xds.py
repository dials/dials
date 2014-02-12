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
  import os
  from libtbx.phil import command_line
  from dials.util.command_line import Importer
  from dxtbx.serialize import xds

  importer = Importer(args)
  experiments = importer.experiments
  reflections = importer.reflections
  args = importer.unhandled_arguments

  if experiments is not None:
    assert len(experiments) > 0

    for i in range(len(experiments)):
      suffix = ""
      if len(experiments) > 1:
        suffix = "_%i" %(i+1)

      sub_dir = "xds%s" %suffix
      if not os.path.isdir(sub_dir):
        os.makedirs(sub_dir)
      # XXX imageset is getting the experimental geometry from the image files
      # rather than the input experiments.json file
      imageset = experiments[i].imageset
      imageset.set_detector(experiments[i].detector)
      imageset.set_beam(experiments[i].beam)
      imageset.set_goniometer(experiments[i].goniometer)
      imageset.set_scan(experiments[i].scan)
      crystal_model = experiments[i].crystal
      crystal_model = crystal_model.change_basis(
        crystal_model.get_space_group().info()\
          .change_of_basis_op_to_reference_setting())
      A = crystal_model.get_A()
      A_inv = A.inverse()
      real_space_a = A_inv.elems[:3]
      real_space_b = A_inv.elems[3:6]
      real_space_c = A_inv.elems[6:9]
      to_xds = xds.to_xds(imageset)
      with open(os.path.join(sub_dir, 'XDS.INP'), 'wb') as f:
        to_xds.XDS_INP(out=f, job_card="XYCORR INIT DEFPIX INTEGRATE CORRECT")
      with open(os.path.join(sub_dir, 'XPARM.XDS'), 'wb') as f:
        to_xds.xparm_xds(
          real_space_a, real_space_b, real_space_c,
          crystal_model.get_space_group().type().number(),
          out=f)

      if reflections is not None and len(reflections) > 0:
        ref_cryst = reflections.select(reflections['id'] == i)
        export_spot_xds(ref_cryst, os.path.join(sub_dir, 'SPOT.XDS'))

  else:
    export_spot_xds(reflections, 'SPOT.XDS')

def export_spot_xds(reflections, filename):
  from iotbx.xds import spot_xds
  if reflections is not None and len(reflections) > 0:
    centroids = reflections['xyzobs.px.value']
    intensities = reflections['intensity.raw.value']
    miller_indices = None
    if 'miller_index' in reflections:
      miller_indices = reflections['miller_index']
      miller_indices = miller_indices.select(miller_indices != (0, 0, 0))
      if len(miller_indices) == 0:
        miller_indices = None
    xds_writer = spot_xds.writer(centroids=centroids,
                                 intensities=intensities,
                                 miller_indices=miller_indices)
    xds_writer.write_file(filename=filename)


if __name__ == '__main__':
  import sys
  run(sys.argv[1:])

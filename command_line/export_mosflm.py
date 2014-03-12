#!/usr/bin/env python
#
# export_mosflm.py
#
#  Copyright (C) 2014 Diamond Light Source
#
#  Author: Richard Gildea
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import division


def run(args):
  import math
  import os
  from scitbx import matrix
  from cctbx.crystal.crystal_model import crystal_model
  from dials.util.command_line import Importer

  importer = Importer(args)
  experiments = importer.experiments
  reflections = importer.reflections
  if reflections is not None:
    assert len(reflections) == 1
    reflections = reflections[0]
  args = importer.unhandled_arguments

  assert len(experiments) > 0

  for i in range(len(experiments)):
    suffix = ""
    if len(experiments) > 1:
      suffix = "_%i" %(i+1)

    sub_dir = "mosflm%s" %suffix
    if not os.path.isdir(sub_dir):
      os.makedirs(sub_dir)
    detector = experiments[i].detector
    beam = experiments[i].beam
    scan = experiments[i].scan
    goniometer = experiments[i].goniometer

    # XXX imageset is getting the experimental geometry from the image files
    # rather than the input experiments.json file
    imageset = experiments[i].imageset

    from rstbx.cftbx.coordinate_frame_helpers import align_reference_frame
    R_to_mosflm = align_reference_frame(
      beam.get_s0(), (1.0, 0.0, 0.0),
      goniometer.get_rotation_axis(), (0.0, 0.0, 1.0))
    #print R_to_mosflm

    cryst = experiments[i].crystal
    cryst = cryst.change_basis(
      cryst.get_space_group().info()\
        .change_of_basis_op_to_reference_setting())
    A = cryst.get_A()
    A_inv = A.inverse()

    real_space_a = R_to_mosflm * A_inv.elems[:3]
    real_space_b = R_to_mosflm * A_inv.elems[3:6]
    real_space_c = R_to_mosflm * A_inv.elems[6:9]

    cryst_mosflm = crystal_model(
      real_space_a, real_space_b, real_space_c,
      space_group=cryst.get_space_group(),
      mosaicity=cryst.get_mosaicity())
    #print cryst_mosflm
    A_mosflm = cryst_mosflm.get_A()
    A_inv_mosflm = A_mosflm.inverse()
    a = matrix.col(A_inv_mosflm[0:3])
    b = matrix.col(A_inv_mosflm[3:6])
    c = matrix.col(A_inv_mosflm[6:9])
    astar = matrix.col(A_mosflm[0:3])
    bstar = matrix.col(A_mosflm[3:6])
    cstar = matrix.col(A_mosflm[6:9])

    a, b, c, alpha, beta, gamma = cryst_mosflm.get_unit_cell().parameters()
    sin_alpha = math.sin(alpha)
    sin_beta = math.sin(beta)
    sin_gamma = math.sin(gamma)
    cos_alpha = math.cos(alpha)
    cos_beta = math.cos(beta)
    cos_gamma = math.cos(gamma)

    # formulae according Rupp A.3 A-45, pg 747.
    V = cryst_mosflm.get_unit_cell().volume()
    V_star = 1/V
    astar = b * c * sin_alpha * V_star
    bstar = a * c * sin_beta * V_star
    cstar = a * b * sin_gamma * V_star
    sin_alpha_star = V / (a * b * c * sin_beta * sin_gamma)
    cos_alpha_star = (cos_beta * cos_gamma - cos_alpha)/ (sin_beta * sin_gamma)
    sin_beta_star = V / (a * b * c * sin_alpha * sin_gamma)
    cos_beta_star = (cos_alpha * cos_gamma - cos_beta)/ (sin_alpha * sin_gamma)
    sin_gamma_star = V / (a * b * c * sin_alpha * sin_beta)
    cos_gamma_star = (cos_alpha * cos_beta - cos_gamma)/ (sin_alpha * sin_beta)

    w = beam.get_wavelength()
    # Mosflm B matrix according to David's notes. This is the same as
    # defined by Busing & Levy apart from the extra factors of 1/w and w.
    B = matrix.sqr((
      astar/w, bstar/w * cos_gamma_star,  cstar/w * cos_beta_star,
            0, bstar/w * sin_gamma_star,     -cstar/w * sin_alpha,
            0,                        0,                      w/c))
    U_mosflm = A_mosflm * B.inverse()

    with open(os.path.join(sub_dir, "index.mat"), "wb") as f:
      print >> f, format_mosflm_mat(w*A_mosflm, U_mosflm, cryst.get_unit_cell())

    directory, template = os.path.split(imageset.get_template())
    symmetry = cryst_mosflm.get_space_group().type().number()
    beam_centre = tuple(reversed(detector[0].get_beam_centre(beam.get_s0())))
    distance = detector[0].get_distance()

    with open(os.path.join(sub_dir, "mosflm.in"), "wb") as f:
      print >> f, write_mosflm_input(directory=directory,
                                     template=template,
                                     symmetry=symmetry,
                                     beam_centre=beam_centre,
                                     distance=distance,
                                     mat_file="index.mat")

  return


def format_mosflm_mat(A, U, unit_cell, missets=(0,0,0)):
  lines = []
  uc_params = unit_cell.parameters()
  for i in range(3):
    lines.append(("%12.8f" * 3) %A.elems[i*3:3*(i+1)])
  lines.append(("%12.3f" * 3) %missets)
  for i in range(3):
    lines.append("%12.8f"*3 %U.elems[i*3:3*(i+1)])
  lines.append(("%12.4f" * 6) %uc_params)
  lines.append(("%12.3f" * 3) %missets)
  return "\n".join(lines)


def write_mosflm_input(directory=None, template=None,
                       symmetry=None,
                       beam_centre=None, distance=None,
                       mat_file=None):
  lines = []
  if directory is not None:
    lines.append("DIRECTORY %s" %directory)
  if template is not None:
    lines.append("TEMPLATE %s" %template)
  if symmetry is not None:
    lines.append("SYMMETRY %s" %symmetry)
  if beam_centre is not None:
    lines.append("BEAM %.3f %.3f" %beam_centre)
  if distance is not None:
    lines.append("DISTANCE %.4f" %distance)
  if mat_file is not None:
    lines.append("MATRIX %s" %mat_file)
  return "\n".join(lines)



if __name__ == '__main__':
  import sys
  run(sys.argv[1:])

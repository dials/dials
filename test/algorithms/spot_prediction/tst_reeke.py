#
#  Copyright (C) (2013) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman.
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#

from __future__ import absolute_import, division
# A class for producing efficient looping limits for reflection
# prediction based on the Reeke algorithm (see Mosflm).

from scitbx import matrix
import scitbx.math
from math import pi
from dials.algorithms.spot_prediction.reeke import reeke_model

def visualize_with_rgl(reeke_model, rscript="reeke_vis.R", dat="reeke_hkl.dat"):
  """Write an R script and an associated data file
  for visualisation of generated indices between phi_beg and phi_end,
  using R and the rgl add-on package."""

  # Sorry, this is ugly. I don't know matplotlib yet.

  # write R script

  with open(rscript, "w") as f:
    f.write("# Run this from within R using\n" + \
            "# install.packages('rgl')\n" + \
            "# source('%s')\n\n" % rscript)
    f.write("library(rgl)\n")
    f.write("p_ax <- c(%.9f, %.9f, %.9f)\n" % reeke_model._rlv_beg[0].elems)
    f.write("q_ax <- c(%.9f, %.9f, %.9f)\n" % reeke_model._rlv_beg[1].elems)
    f.write("r_ax <- c(%.9f, %.9f, %.9f)\n" % reeke_model._rlv_beg[2].elems)
    f.write("source <- c(%.9f, %.9f, %.9f)\n" % reeke_model._source.elems)
    f.write("rot_ax <- c(%.9f, %.9f, %.9f)\n" % reeke_model._axis.elems)
    f.write("phi_range <- c(%.9f, %.9f)\n" % reeke_model._phi_range)
    f.write("half_osc <- matrix(data = c(" + \
          "%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f" % \
          reeke_model._r_half_osc.elems + \
          "), nrow=3, byrow=T)\n")
    f.write("dstarmax <- %.9f\n" % reeke_model._dstarmax)
    f.write("perm <- solve(matrix(data = c(" + \
          "%d,%d,%d,%d,%d,%d,%d,%d,%d" % reeke_model._permutation.elems + \
          "), nrow=3, byrow=T))\n")
    f.write("\n# draw the Ewald and limiting spheres\n" + \
          "open3d()\n" + \
          "spheres3d(source,radius=sqrt(sum(source*source)),color='#CCCCFF'," + \
          "alpha=0.3)\n" + \
          "spheres3d(c(0,0,0),radius=dstarmax," + \
          "color='red',alpha=0.1)\n" + \
          "\n# draw the source vector and rotation axis\n" + \
          "lines3d(rbind(c(0,0,0),source))\n" + \
          "lines3d(rbind(c(0,0,0),rot_ax))\n" + \
          "\n# draw the reciprocal lattice axes at ten times their " + \
          "length\n" + \
          "lines3d(rbind(c(0,0,0),10*p_ax),col='red')\n" + \
          "lines3d(rbind(c(0,0,0),10*q_ax),col='green')\n" + \
          "lines3d(rbind(c(0,0,0),10*r_ax),col='blue')\n"
          "sourcemag <- sqrt(sum(source*source))\n" + \
          "sourceunit <- source / sourcemag\n" + \
          "\n# two unit vectors orthogonal to source\n" + \
          "sourceunitx1 <- c((-sourceunit[3]),(0),(sourceunit[1]))\n" + \
          "sourceunitx2 <- c((sourceunit[2]*sourceunitx1[3] - " + \
          "sourceunit[3]*sourceunitx1[2]),\n" + \
          "   (sourceunit[1]*sourceunitx1[3] - sourceunit[3]*sourceunitx1[1]),\n" + \
          "   (sourceunit[1]*sourceunitx1[2] - sourceunit[2]*sourceunitx1[1]))\n" + \
          "sin_theta <- dstarmax/(2*sourcemag)\n" + \
          "sin_2theta <- sin(2*asin(sin_theta))\n" + \
          "\n# distance to the centre of the circle of" + \
          "intersection, along source\n" + \
          "e <- 2 * sqrt(sum(source*source)) * sin_theta ^2\n" + \
          "\n# radius of the circle of intersection\n" + \
          "R <- sourcemag * sin_2theta\n" + \
          "\n# make points around the circle\n" + \
          "tau <- seq(from=0,to=2*pi,by=0.01)\n" + \
          "circ <- t(sapply(tau, function(x){\n" + \
          "  e * sourceunit + R*sin(x) * sourceunitx1 + R*cos(x) * sourceunitx2}))\n" + \
          "\n# draw the circle\n" + \
          "lines3d(circ)\n" + \
          "\n# load the generated indices\n"
          "pts <- read.csv('./%s')\n" % dat)
    f.write("\n# convert h, k, l to reciprocal space coordinates\n" + \
          "conv <- function(h) {p <- perm %*% h\n" + \
          "    return(p[1]*p_ax + p[2]*q_ax + p[3]*r_ax)}\n" + \
          "pts <- t(apply(pts, MARGIN = 1, FUN = conv))\n" + \
          "\n# draw the generated indices\n" + \
          "points3d(pts, col='blue')\n" + \
          "\n")

  # write data file

  indices = reeke_model.generate_indices()

  with open(dat, "w") as f:
    f.write("h, k, l\n")
    for hkl in indices:
      f.write("%d, %d, %d\n" % hkl)

  print "Generated indices were written to %s" % dat
  print "An R script for visualising these was written to %s," % rscript
  print "which can be run from the R prompt with:"
  print "source('%s')" % rscript

  return

def reeke_model_for_use_case(phi_beg, phi_end, margin):
  """Construct a reeke_model for the geometry of the Use Case Thaumatin
  dataset, taken from the XDS XPARM. The values are hard-
  coded here so that this module does not rely on the location of that
  file."""

  axis = matrix.col([0.0, 1.0, 0.0])

  # original (unrotated) setting
  ub = matrix.sqr([-0.0133393674072, -0.00541609051856, -0.00367748834997,
                  0.00989309470346, 0.000574825936669, -0.0054505379664,
                  0.00475395109417, -0.0163935257377, 0.00102384915696])
  r_beg = matrix.sqr(scitbx.math.r3_rotation_axis_and_angle_as_matrix(
      axis = self._axis, angle = phi_beg, deg = True))
  r_osc = matrix.sqr(
      scitbx.math.r3_rotation_axis_and_angle_as_matrix(
      axis = self._axis, angle = (phi_end - phi_beg), deg=True))

  ub_beg = r_beg * ub
  ub_end = self._r_osc * ub_mid
  s0 = matrix.col([0.00237878589035, 1.55544539299e-16, -1.09015329696])
  dmin = 1.20117776325

  return reeke_model(ub_beg, ub_end, axis, s0, dmin, margin)

def regression_test():
  """Perform a regression test by comparing to indices generating
  by the brute force method used in the Use Case."""

  from rstbx.diffraction import rotation_angles
  from rstbx.diffraction import full_sphere_indices
  from cctbx.sgtbx import space_group, space_group_symbols
  from cctbx.uctbx import unit_cell

  # cubic, 50A cell, 1A radiation, 1 deg osciillation, everything ideal

  a = 50.0

  ub_beg = matrix.sqr((1.0 / a, 0.0, 0.0,
                       0.0, 1.0 / a, 0.0,
                       0.0, 0.0, 1.0 / a))

  axis = matrix.col((0, 1, 0))

  r_osc = matrix.sqr(
          scitbx.math.r3_rotation_axis_and_angle_as_matrix(
          axis = axis, angle = 1.0, deg=True))

  ub_end = r_osc * ub_beg

  uc = unit_cell((a, a, a, 90, 90, 90))
  sg = space_group(space_group_symbols('P23').hall())

  s0 = matrix.col((0, 0, 1))

  wavelength = 1.0
  dmin = 1.5

  indices = full_sphere_indices(
      unit_cell = uc, resolution_limit = dmin, space_group = sg)

  ra = rotation_angles(dmin, ub_beg, wavelength, axis)

  obs_indices, obs_angles = ra.observed_indices_and_angles_from_angle_range(
      phi_start_rad = 0.0 * pi / 180.0,
      phi_end_rad = 1.0 * pi / 180.0,
      indices = indices)

  r = reeke_model(ub_beg, ub_end, axis, s0, dmin, 1.0)
  reeke_indices = r.generate_indices()
  #r.visualize_with_rgl()

  for oi in obs_indices:
    assert(tuple(map(int, oi)) in reeke_indices)

  #TODO Tests for an oblique cell


if __name__ == '__main__':

  import sys

  if len(sys.argv) == 1:
    regression_test()

  elif len(sys.argv) < 3:
    from libtbx.utils import Sorry
    raise Sorry("Expecting either 3 or 4 arguments: path/to/xparm.xds start_phi end_phi margin=3")

  else:

    # take an xparm.xds, phi_beg, phi_end and margin from the command arguments.
    from rstbx.cftbx.coordinate_frame_converter import \
        coordinate_frame_converter
    cfc = coordinate_frame_converter(sys.argv[1])
    phi_beg, phi_end = float(sys.argv[2]), float(sys.argv[3])
    margin = int(sys.argv[4]) if len(sys.argv) == 5 else 3

    # test run for development/debugging.
    u, b = cfc.get_u_b()
    ub = matrix.sqr(u * b)
    axis = matrix.col(cfc.get('rotation_axis'))
    rub_beg = axis.axis_and_angle_as_r3_rotation_matrix(phi_beg) * ub
    rub_end = axis.axis_and_angle_as_r3_rotation_matrix(phi_end) * ub
    wavelength = cfc.get('wavelength')
    sample_to_source_vec = matrix.col(cfc.get_c('sample_to_source').normalize())
    s0 = (- 1.0 / wavelength) * sample_to_source_vec
    dmin = 1.20117776325

    r = reeke_model(rub_beg, rub_end, axis, s0, dmin, margin)

    indices = r.generate_indices()

    for hkl in indices:
      print "%4d %4d %4d" % hkl

#
#  Copyright (C) (2013) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman.
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#

"""Auxiliary functions for the refinement package"""

from __future__ import division
from math import sin, cos, acos
from scitbx import matrix
from scitbx.array_family import flex #import dependency
from dials_refinement_helpers_ext import dR_from_axis_and_angle as dR_cpp
from dials_refinement_helpers_ext import CrystalOrientationCompose as xloc_cpp
import random

def ordinal_number(array_index=None, cardinal_number=None):
  '''Return a string representing the ordinal number for the input integer. One
  of array_index or cardinal_number must be set, depending on whether the
  input is from a 0-based or 1-based sequence.

  Based on Thad Guidry's post at
  https://groups.google.com/forum/#!topic/openrefine/G7_PSdUeno0'''
  if [array_index, cardinal_number].count(None) != 1:
    raise ValueError("One of array_index or cardinal_number should be set")
  if array_index is not None:
    i = int(array_index) + 1
  if cardinal_number is not None:
    i = int(cardinal_number)
  return str(i) + {1: 'st', 2: 'nd', 3: 'rd'}.get(4 if 10 <= i % 100 < 20 else i % 10, "th")

class CrystalOrientationCompose(xloc_cpp):
  '''Wrapper for the C++ CrystalOrientationCompose class wiht accessors that
  return matrix.sqr values.'''

  def U(self):
    return matrix.sqr(super(CrystalOrientationCompose, self).U())

  def dU_dphi1(self):
    return matrix.sqr(super(CrystalOrientationCompose, self).dU_dphi1())

  def dU_dphi2(self):
    return matrix.sqr(super(CrystalOrientationCompose, self).dU_dphi2())

  def dU_dphi3(self):
    return matrix.sqr(super(CrystalOrientationCompose, self).dU_dphi3())

def dR_from_axis_and_angle(axis, angle, deg=False):
  """Wrapper for C++ version of dR_from_axis_and_angle returning a matrix.sqr"""
  return matrix.sqr(dR_cpp(axis, angle, deg))

def dR_from_axis_and_angle_py(axis, angle, deg=False):
  """return the first derivative of a rotation matrix specified by its
  axis and angle"""

  # NB it is inefficient to do this separately from the calculation of
  # the rotation matrix itself, but it seems the Python interface to
  # scitbx does not have a suitable function. It might perhaps be
  # useful to write one, which could come straight from David Thomas'
  # RTMATS (present in Mosflm and MADNES).

  # NB RTMATS does calculation for a clockwise rotation of a vector
  # whereas axis_and_angle_as_r3_rotation_matrix does anticlockwise
  # rotation. Therefore flip the axis in here compared with
  # RTMATS in order to match the axis_and_angle_as_r3_rotation_matrix
  # convention

  # See also axis_and_angle_as_r3_derivative_wrt_angle, which does the same
  # as this function, but this function is faster.

  assert axis.n in ((3,1), (1,3))
  if (deg): angle *= pi/180
  axis = -1. * axis.normalize()
  ca, sa  = cos(angle), sin(angle)

  return(matrix.sqr((sa * axis[0] * axis[0] - sa ,
                  sa * axis[0] * axis[1] + ca * axis[2],
                  sa * axis[0] * axis[2] - ca * axis[1],
                  sa * axis[1] * axis[0] - ca * axis[2],
                  sa * axis[1] * axis[1] - sa,
                  sa * axis[1] * axis[2] + ca * axis[0],
                  sa * axis[2] * axis[0] + ca * axis[1],
                  sa * axis[2] * axis[1] - ca * axis[0],
                  sa * axis[2] * axis[2] - sa)))

def skew_symm(v):
  '''Make matrix [v]_x from v. Essentially multiply vector by SO(3) basis
  set Lx, Ly, Lz. Equation (2) from Gallego & Yezzi paper.

  NB a C++ version exists in gallego_yezzi.h.'''
  import scitbx.matrix

  L1 = scitbx.matrix.sqr((0, 0, 0, 0, 0, -1, 0, 1, 0))
  L2 = scitbx.matrix.sqr((0, 0, 1, 0, 0, 0, -1, 0, 0))
  L3 = scitbx.matrix.sqr((0, -1, 0, 1, 0, 0, 0, 0, 0))

  v1, v2, v3 = v.elems

  return v1 * L1 + v2 * L2 + v3 * L3

def dRq_de(theta, e, q):
  '''Calculate derivative of rotated vector r = R*q with respect to the elements
  of the rotation axis e, where the angle of rotation is theta.

  Implementation of Equation (8) from Gallego & Yezzi.

  NB a C++ version exists in gallego_yezzi.h.'''

  from scitbx import matrix

  # ensure e is unit vector
  e = e.normalize()

  # rotation matrix
  R = e.axis_and_angle_as_r3_rotation_matrix(theta, deg=False)

  # rotation vector v
  v = theta * e

  qx = skew_symm(q)
  vx = skew_symm(v)
  vvt = v*v.transpose()
  Rt = R.transpose()
  I3 = matrix.identity(3)

  return (-1./theta) * R * qx * (vvt + (Rt - I3) * vx)

def random_param_shift(vals, sigmas):
  """Add a random (normal) shift to a parameter set, for testing"""

  assert len(vals) == len(sigmas)
  shifts = [random.gauss(0, sd) for sd in sigmas]
  newvals = [(x + y) for x, y in zip(vals, shifts)]

  return newvals

def get_fd_gradients(mp, deltas, multi_state_elt=None):
  """Calculate centered finite difference gradients for each of the
  parameters of the model parameterisation mp.

  "deltas" must be a sequence of the same length as the parameter list, and
  contains the step size for the difference calculations for each parameter.

  "multi_state_elt" selects a particular state for use when mp is a multi-
  state parameterisation.
  """

  p_vals = mp.get_param_vals()
  assert len(deltas) == len(p_vals)
  fd_grad = []

  for i in range(len(deltas)):

    val = p_vals[i]

    p_vals[i] -= deltas[i] / 2.
    mp.set_param_vals(p_vals)
    if multi_state_elt is None:
      rev_state = mp.get_state()
    else:
      rev_state = mp.get_state(multi_state_elt=multi_state_elt)

    p_vals[i] += deltas[i]
    mp.set_param_vals(p_vals)
    if multi_state_elt is None:
      fwd_state = mp.get_state()
    else:
      fwd_state = mp.get_state(multi_state_elt=multi_state_elt)

    fd_grad.append((fwd_state - rev_state) / deltas[i])

    p_vals[i] = val

  # return to the initial state
  mp.set_param_vals(p_vals)

  return fd_grad

def get_panel_groups_at_depth(group, depth=0):
  """Return a list of the panel groups at a certain depth below the node group"""
  assert depth >= 0
  if depth == 0:
    return [group]
  else:
    return [p for gp in group.children() for p in get_panel_groups_at_depth(gp, depth-1)]

def get_panel_ids_at_root(panel_list, group):
  """Get the sequential panel IDs for a set of panels belonging to a group"""
  try:
    return [p for gp in group.children() for p in get_panel_ids_at_root(panel_list, gp)]
  except AttributeError: # we got down to Panels
    return [panel_list.index(group)]

def corrgram(corrmat, labels):
  """Create a correlation matrix plot or 'corrgram' for the provided
  correlation matrix and row/column labels. Inspired by R's corrplot and
  https://github.com/louridas/corrplot/blob/master/corrplot.py"""

  try: # is corrmat a scitbx matrix?
    corrmat = corrmat.as_flex_double_matrix()
  except AttributeError: # assume it is already a flex double matrix
    pass
  assert corrmat.is_square_matrix()

  nr = corrmat.all()[0]
  assert nr == len(labels)

  from math import pi, sqrt
  try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
  except ImportError as e:
    msg = "matplotlib modules not available " + str(e)
    info(msg)
    return None

  plt.figure(1)
  ax = plt.subplot(1, 1, 1, aspect='equal')
  clrmap = cm.get_cmap('bwr')

  for x in xrange(nr):
    for y in xrange(nr):
      d = corrmat[x, y]
      d_abs = abs(d)
      circ = plt.Circle((x, y),radius=0.9*sqrt(d_abs)/2)
      circ.set_edgecolor('white')
      # put data into range [0,1] and invert so that 1 == blue and 0 == red
      facecolor = 1 - (0.5*d + 0.5)
      circ.set_facecolor(clrmap(facecolor))
      ax.add_artist(circ)
  ax.set_xlim(-0.5, nr-0.5)
  ax.set_ylim(-0.5, nr-0.5)

  ax.xaxis.tick_top()
  xtickslocs = range(len(labels))
  ax.set_xticks(xtickslocs)
  ax.set_xticklabels(labels, rotation=30, fontsize='small', ha='left')

  ax.invert_yaxis()
  ytickslocs = range(len(labels))
  ax.set_yticks(ytickslocs)
  ax.set_yticklabels(labels, fontsize='small')

  xtickslocs = [e + 0.5 for e in range(len(labels))]
  ax.set_xticks(xtickslocs, minor=True)
  ytickslocs = [e + 0.5 for e in range(len(labels))]
  ax.set_yticks(ytickslocs, minor=True)
  plt.grid(color='0.8', which='minor', linestyle='-')

  # suppress major tick marks
  ax.tick_params(which='major', width=0)

  # need this otherwise text gets clipped
  plt.tight_layout()

  # FIXME should this also have a colorbar as legend?
  return plt

def string_sel(l, full_names, prefix=""):
  '''Provide flexible matching between a list of input strings, l,
  consisting either of indices or partial names, and a list of full names,
  with an optional shared prefix. The input list l may arrive from PHIL
  conversion of the strings type. In that case, comma-separated values will
  require splitting, and bracket characters will be removed. The values in
  the processed list l should consist of integers or partial names. Integers
  will be treated as 0-based indices and partial names will be matched to
  as many full names as possible. The created selection is returned as a
  boolean list.'''

  sel = [False] * len(full_names)
  full_names = [prefix + s for s in full_names]

  # expand elements of the list that are comma separated strings and remove
  # braces/brackets
  l = [s.strip('(){}[]') for e in l for s in str(e).split(',')]
  l = [e for e in l if e is not '']
  for e in l:
    try:
      i = int(e)
      sel[i] = True
      continue
    except ValueError:
      pass
    except IndexError:
      pass
    sel = [True if e in name else s for (name, s) in \
      zip(full_names, sel)]

  return sel

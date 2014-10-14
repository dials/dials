from __future__ import division

import math
from cctbx.array_family import flex
import iotbx.phil
from cctbx import crystal, miller
from scitbx import matrix

master_phil_scope = iotbx.phil.parse(
"""
hkl = None
  .type = ints(size=3)
  .multiple=True
hkl_limit = None
  .type = int(value_min=1)
expand_to_p1 = True
  .type = bool
  .help = "Expand the given miller indices to symmetry equivalent reflections"
eliminate_sys_absent = True
  .type = bool
  .help = "Eliminate systematically absent reflections"
plane_normal = 0,0,1
  .type = ints(size=3)
save_coordinates = True
  .type = bool
plot = True
  .type = bool
""")


def stereographic_projection(crystal_model, miller_indices,
                             plane_normal=(0,0,1),
                             rotation_matrix=None):
  # http://dx.doi.org/10.1107/S0021889868005029
  # J. Appl. Cryst. (1968). 1, 68-70
  # The construction of stereographic projections by computer
  # G. K. Stokes, S. R. Keown and D. J. Dyson

  A = crystal_model.get_A()
  A_inv = A.inverse()

  G = A_inv * A_inv.transpose()
  G_star = A.transpose() * A
  from libtbx.test_utils import approx_equal
  assert approx_equal(G_star, G.inverse())
  h0 = (G * matrix.col(plane_normal)).normalize()
  h1 = matrix.col((1,0,0)).cross((G_star * h0).normalize())
  h2 = (G_star * h1).cross(G_star * h0).normalize()

  if rotation_matrix is not None:
    h0 = rotation_matrix * h0
    h1 = rotation_matrix * h1
    h2 = rotation_matrix * h2

  h0_t = h0.transpose()
  h1_t = h1.transpose()
  h2_t = h2.transpose()

  projections = flex.vec2_double()

  for hkl in miller_indices:
    hi = matrix.col(hkl)
    hi_t = hi.transpose()
    for sign in (+1, -1):
      hi *= sign
      hi_t *= sign
      cos_theta = (hi_t * G_star * h0)[0]/(
        math.sqrt((hi_t * G_star * hi)[0] * (h0_t * G_star * h0)[0]))
      if sign == -1:
        hi *= sign
        hi_t *= sign
      if cos_theta >= 0:
        break
      else:
        continue
    cos_alpha = (hi_t * G_star * h1)[0]/(
      math.sqrt((hi_t * G_star * hi)[0] * (h1_t * G_star * h1)[0]))
    N = (hi_t * G_star * h2)[0]
    if abs(cos_theta) > 1:
      cos_theta = math.copysign(1, cos_theta) # XXX why?
    theta = math.acos(cos_theta)
    if theta == 0:
      projections.append((0,0)) # XXX
      continue
    r = math.tan(theta/2)
    cos_phi = cos_alpha / math.sin(theta)
    if abs(cos_phi) > 1:
      cos_phi = math.copysign(1, cos_phi) # XXX why?
    x = r * cos_phi
    y = math.copysign(abs(r * math.sin(math.acos(cos_phi))), N)
    projections.append((x,y))
  return projections


def gcd_list(l):
  # greatest common divisor for a list of numbers
  from scitbx.math import gcd_int_simple as gcd
  result = l[0]
  for i in range(1, len(l)):
    result = gcd(result, l[i])
  return result


def run(args):
  #print args

  from dials.util.options import OptionParser
  from dials.util.options import flatten_experiments

  parser = OptionParser(
    phil=master_phil_scope,
    read_experiments=True,
    check_format=False)

  params, options = parser.parse_args(show_diff_phil=True)
  experiments = flatten_experiments(params.input.experiments)
  assert len(experiments) > 0

  if len(params.hkl) == 0 and params.hkl_limit is None:
    from libtbx.utils import Sorry
    raise Sorry("Please provide hkl or hkl_limit parameters.")

  if params.hkl is not None and len(params.hkl):
    miller_indices = flex.miller_index(params.hkl)
  elif params.hkl_limit is not None:
    limit = params.hkl_limit
    miller_indices = flex.miller_index()
    for h in range(-limit, limit+1):
      for k in range(-limit, limit+1):
        for l in range(-limit, limit+1):
          if (h,k,l) == (0,0,0): continue
          miller_indices.append((h,k,l))

  crystals = experiments.crystals()

  symmetry = crystal.symmetry(
    unit_cell=crystals[0].get_unit_cell(),
    space_group=crystals[0].get_space_group())
  miller_set = miller.set(symmetry, miller_indices)
  d_spacings = miller_set.d_spacings()
  if params.eliminate_sys_absent:
    d_spacings = d_spacings.eliminate_sys_absent()
  if params.expand_to_p1:
    d_spacings = d_spacings.as_non_anomalous_array().expand_to_p1()
    d_spacings = d_spacings.generate_bijvoet_mates()
  miller_indices = d_spacings.indices()

  # find the greatest common factor (divisor) between miller indices
  miller_indices_unique = flex.miller_index()
  for hkl in miller_indices:
    gcd = gcd_list(hkl)
    if gcd > 1:
      miller_indices_unique.append(tuple(int(h/gcd) for h in hkl))
    elif gcd < 1:
      pass
    else:
      miller_indices_unique.append(hkl)
  miller_indices = miller_indices_unique
  miller_indices = flex.miller_index(list(set(miller_indices)))

  plane_normal = matrix.col(params.plane_normal).normalize()

  ref_crystal = crystals[0]
  projections_ref = stereographic_projection(
    ref_crystal, miller_indices, plane_normal=plane_normal.elems)

  projections_all = [projections_ref]

  if len(crystals) > 0:
    from dials.algorithms.indexing.compare_orientation_matrices import \
        difference_rotation_matrix_and_euler_angles

    for cryst in crystals[1:]:
      R_ij, euler_angles, cb_op = difference_rotation_matrix_and_euler_angles(
        ref_crystal, cryst)
      projections = stereographic_projection(
        cryst, miller_indices, plane_normal=plane_normal.elems,
        rotation_matrix=R_ij)
      projections_all.append(projections)

  if params.plot:
    from matplotlib import pyplot
    from matplotlib import pylab

    cir = pylab.Circle((0,0), radius=1.0, fill=False, color='0.75')
    pylab.gca().add_patch(cir)

    for projections in projections_all:
      x, y = projections.parts()
      pyplot.scatter(x.as_numpy_array(), y.as_numpy_array(), c='b', s=2)
      #for hkl, proj in zip(miller_indices, projections):
        #pyplot.text(proj[0], proj[1], str(hkl), fontsize=5)
    pyplot.axes().set_aspect('equal')
    pyplot.xlim(-1.1,1.1)
    pyplot.ylim(-1.1,1.1)
    pyplot.show()

  if params.save_coordinates:
    with open('projections.txt', 'wb') as f:
      for projections in projections_all:
        for hkl, proj in zip(miller_indices, projections):
          print >> f, "%i %i %i" %hkl,
          print >> f, "%f %f" %proj



if __name__ == '__main__':
  import sys
  run(sys.argv[1:])

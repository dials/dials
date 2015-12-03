"""
Figure out correct gradient expressions required for crystal unit cell
restraints

"""

# Python and cctbx imports
from __future__ import division
import sys
from math import pi, sqrt
from scitbx import matrix
from libtbx.phil import parse
from libtbx.test_utils import approx_equal
from scitbx.array_family import flex

# Get modules to build models and minimiser using PHIL
from dials.test.algorithms.refinement import setup_geometry
from dials.test.algorithms.refinement import setup_minimiser

from dials.algorithms.refinement.parameterisation.crystal_parameters import \
    CrystalOrientationParameterisation, CrystalUnitCellParameterisation

# Symmetry constrained parameterisation for the unit cell
from cctbx.uctbx import unit_cell
from rstbx.symmetry.constraints.parameter_reduction import \
    symmetrize_reduce_enlarge

DEG2RAD = pi/180.0
RAD2DEG = 180.0/pi

#args = sys.argv[1:]
master_phil = parse("""
    include scope dials.test.algorithms.refinement.geometry_phil
    include scope dials.test.algorithms.refinement.minimiser_phil
    """, process_includes=True)

# make cell more oblique
args=["a.direction.close_to.sd=5","b.direction.close_to.sd=5","c.direction.close_to.sd=5"]
models = setup_geometry.Extract(master_phil, cmdline_args = args)
crystal = models.crystal

print crystal

# derive finite difference gradients of various quantities wrt each param
def check_fd_gradients(parameterisation):

  mp = parameterisation
  p_vals = mp.get_param_vals()
  deltas = [1.e-7 for p in p_vals]
  assert len(deltas) == len(p_vals)
  fd_grad = []

  for i in range(len(deltas)):

    val = p_vals[i]

    p_vals[i] -= deltas[i] / 2.
    mp.set_param_vals(p_vals)
    rev_uc = mp.get_model().get_unit_cell().parameters()
    rev_vec = mp.get_model().get_real_space_vectors()
    rev_B = mp.get_model().get_B()
    rev_O = rev_B.transpose().inverse()

    p_vals[i] += deltas[i]
    mp.set_param_vals(p_vals)
    fwd_uc = mp.get_model().get_unit_cell().parameters()
    fwd_vec = mp.get_model().get_real_space_vectors()
    fwd_B = mp.get_model().get_B()
    fwd_O = fwd_B.transpose().inverse()

    fd_uc = [(f - r) / deltas[i] for f,r in zip(fwd_uc, rev_uc)]
    fd_vec = [(f - r) / deltas[i] for f,r in zip(fwd_vec, rev_vec)]
    fd_B = (fwd_B - rev_B) / deltas[i]
    fd_O = (fwd_O - rev_O) / deltas[i]
    fd_grad.append({'da_dp':fd_uc[0],
                    'db_dp':fd_uc[1],
                    'dc_dp':fd_uc[2],
                    'daa_dp':fd_uc[3],
                    'dbb_dp':fd_uc[4],
                    'dcc_dp':fd_uc[5],
                    'davec_dp':fd_vec[0],
                    'dbvec_dp':fd_vec[1],
                    'dcvec_dp':fd_vec[2],
                    'dB_dp':fd_B,
                    'dO_dp':fd_O})

    p_vals[i] = val

  # return to the initial state
  mp.set_param_vals(p_vals)

  return fd_grad

xlo_param = CrystalOrientationParameterisation(crystal)
xluc_param = CrystalUnitCellParameterisation(crystal)

from dials.algorithms.refinement.restraints.restraints import SingleUnitCellTie
uct = SingleUnitCellTie(xluc_param, [None]*6, [None]*6)

from dials.algorithms.refinement.refinement_helpers import \
      AngleDerivativeWrtVectorElts

B = crystal.get_B()
O = (B.transpose()).inverse()
a, b, c, aa, bb, cc = crystal.get_unit_cell().parameters()
aa *= DEG2RAD
bb *= DEG2RAD
cc *= DEG2RAD
avec, bvec, cvec = crystal.get_real_space_vectors()

# calculate d[B^T]/dp
dB_dp = xluc_param.get_ds_dp()
dBT_dp = [dB.transpose() for dB in dB_dp]

# calculate d[O]/dp
dO_dp = [-O * dBT * O for dBT in dBT_dp]

# objects to get derivative of angles wrt vectors
dalpha = AngleDerivativeWrtVectorElts(bvec, cvec)
dbeta = AngleDerivativeWrtVectorElts(avec, cvec)
dgamma = AngleDerivativeWrtVectorElts(avec, bvec)

# get all FD derivatives
fd_grad = check_fd_gradients(xluc_param)

# look at each parameter
for i, dO in enumerate(dO_dp):

  #print "dB_dp analytical"
  #print dB_dp[i]
  #print "dB_dp FD"
  #print fd_grad[i]['dB_dp']
  #print

  print "***** PARAMETER {0} *****".format(i)

  # dB_dp is good.

  print "O MATRIX"
  print "dO_dp analytical"
  print dO
  print "dO_dp FD"
  print fd_grad[i]['dO_dp']
  print

  # extract derivatives of each unit cell vector wrt p
  dav_dp, dbv_dp, dcv_dp = dO.transpose().as_list_of_lists()
  dav_dp = matrix.col(dav_dp)
  dbv_dp = matrix.col(dbv_dp)
  dcv_dp = matrix.col(dcv_dp)

  # check these are correct vs FD
  print "CELL VECTORS"
  diff = dav_dp - fd_grad[i]['davec_dp']
  #print 2 * diff.length() / (dav_dp.length() + fd_grad[i]['davec_dp'].length()) * 100
  print 'davec_dp analytical: {0} {1} {2}'.format(*dav_dp.elems)
  print 'davec_dp finite diff: {0} {1} {2}'.format(*fd_grad[i]['davec_dp'].elems)

  # only the first one seems about right. What about b?
  diff = dbv_dp - fd_grad[i]['dbvec_dp']
  #print 2 * diff.length() / (dbv_dp.length() + fd_grad[i]['dbvec_dp'].length()) * 100
  print 'dbvec_dp analytical: {0} {1} {2}'.format(*dbv_dp.elems)
  print 'dbvec_dp finite diff: {0} {1} {2}'.format(*fd_grad[i]['dbvec_dp'].elems)

  # and c?
  diff = dcv_dp - fd_grad[i]['dcvec_dp']
  #print 2 * diff.length() / (dcv_dp.length() + fd_grad[i]['dcvec_dp'].length()) * 100
  print 'dcvec_dp analytical: {0} {1} {2}'.format(*dcv_dp.elems)
  print 'dcvec_dp finite diff: {0} {1} {2}'.format(*fd_grad[i]['dcvec_dp'].elems)
  print

  #
  print "CELL LENGTHS"
  da_dp = 1./a * avec.dot(dav_dp)
  print "d[a]/dp{2} analytical: {0} FD: {1}".format(da_dp, fd_grad[i]['da_dp'], i)

  db_dp = 1./b * bvec.dot(dbv_dp)
  print "d[b]/dp{2} analytical: {0} FD: {1}".format(db_dp, fd_grad[i]['db_dp'], i)

  dc_dp = 1./c * cvec.dot(dcv_dp)
  print "d[c]/dp{2} analytical: {0} FD: {1}".format(dc_dp, fd_grad[i]['dc_dp'], i)

  # dc_dp appears worse than da_dp or db_dp for the 6th parameter. The vector dcv_dp
  # is 20 times longer than dbv_dp, so clearly the c vector is strongly affected
  # by this parameter, which is g5=b*.c*. However, the direction of this change
  # is almost perpendicular to cvec.

  print
  print "CELL ANGLES"

  #z = bvec.dot(cvec) / (b * c)
  #daa_dp = bvec.dot(cvec) * (db_dp * c + b * dc_dp) - b * c * (dbv_dp.dot(cvec) + bvec.dot(dcv_dp))
  #daa_dp /= (b * b * c * c)
  #daa_dp *= -RAD2DEG / (sqrt(1 - z**2))

  # Here we know the derivatives of the angle alpha with respect to elements
  # of the vectors a and b. We know these expressions are correct because they
  # are tested in tst_angle_derivatives_wrt_vector_elts.py
  dalpha_db = dalpha.derivative_wrt_u()
  dalpha_dc = dalpha.derivative_wrt_v()
  dbeta_da = dbeta.derivative_wrt_u()
  dbeta_dc = dbeta.derivative_wrt_v()
  dgamma_da = dgamma.derivative_wrt_u()
  dgamma_db = dgamma.derivative_wrt_v()

  # why is this wrong?
  daa_dp = RAD2DEG * dbv_dp.dot(dalpha_db) + dcv_dp.dot(dalpha_dc)
  dbb_dp = RAD2DEG * dav_dp.dot(dbeta_da) + dcv_dp.dot(dbeta_dc)
  dcc_dp = RAD2DEG * dav_dp.dot(dgamma_da) + dbv_dp.dot(dgamma_db)

  # because only orthogonal changes are relevant?
  ua = avec.normalize()
  if dav_dp.length() < 1e-10:
    ortho_dav_dp = matrix.col((0, 0, 0))
  else:
    v = avec.cross(dav_dp).normalize()
    u = v.cross(ua).normalize()
    ortho_dav_dp = dav_dp.dot(u) * u

  ub = bvec.normalize()
  if dbv_dp.length() < 1e-10:
    ortho_dbv_dp = matrix.col((0, 0, 0))
  else:
    v = bvec.cross(dbv_dp).normalize()
    u = v.cross(ub).normalize()
    ortho_dbv_dp = dbv_dp.dot(u) * u

  uc = cvec.normalize()
  if dcv_dp.length() < 1e-10:
    ortho_dcv_dp = matrix.col((0, 0, 0))
  else:
    v = cvec.cross(dcv_dp).normalize()
    u = v.cross(uc).normalize()
    ortho_dcv_dp = dcv_dp.dot(u) * u

  daa_dp = RAD2DEG * (ortho_dbv_dp.dot(dalpha_db) + ortho_dcv_dp.dot(dalpha_dc))
  dbb_dp = RAD2DEG * (ortho_dav_dp.dot(dbeta_da) + ortho_dcv_dp.dot(dbeta_dc))
  dcc_dp = RAD2DEG * (ortho_dav_dp.dot(dgamma_da) + ortho_dbv_dp.dot(dgamma_db))

  print "d[alpha]/dp{2} analytical: {0} FD: {1}".format(daa_dp, fd_grad[i]['daa_dp'], i)
  print "d[beta]/dp{2} analytical: {0} FD: {1}".format(dbb_dp, fd_grad[i]['dbb_dp'], i)
  print "d[gamma]/dp{2} analytical: {0} FD: {1}".format(dcc_dp, fd_grad[i]['dcc_dp'], i)

  #print "analytical dcc_dp", dcc_dp * RAD2DEG
  #print "FD dcc_dp", fd_grad[i]['dcc_dp']
  #print

# enter interactive console
#from dials.util.command_line import interactive_console; interactive_console()

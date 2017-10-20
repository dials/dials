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

master_phil = parse("""
    include scope dials.test.algorithms.refinement.geometry_phil
    include scope dials.test.algorithms.refinement.minimiser_phil
    """, process_includes=True)

# make cell more oblique
args=["a.direction.close_to.sd=5","b.direction.close_to.sd=5","c.direction.close_to.sd=5"]
models = setup_geometry.Extract(master_phil, cmdline_args = args)
crystal = models.crystal

# a hexagonal crystal is a good test case for behaviour of oblique cells
do_hexagonal = True
if do_hexagonal:
  import libtbx.load_env, os
  from dxtbx.model.experiment_list import ExperimentListFactory
  try:
    dials_regression = libtbx.env.dist_path('dials_regression')
  except KeyError:
    print 'Skipped: dials_regression not configured'
    exit(0)
  data_dir = os.path.join(dials_regression, "refinement_test_data", "multi_stills")
  experiments_path = os.path.join(data_dir, "combined_experiments.json")
  experiments = ExperimentListFactory.from_json_file(experiments_path,
                check_format=False)
  crystal = experiments[0].crystal

# derive finite difference gradients of various quantities wrt each param
def check_fd_gradients(parameterisation):

  mp = parameterisation
  p_vals = mp.get_param_vals()
  deltas = [1.e-7 for p in p_vals]
  assert len(deltas) == len(p_vals)
  fd_grad = []

  # get matrix to unset rotations of unit cell vectors
  Ut = matrix.sqr(mp.get_model().get_U()).transpose()

  for i in range(len(deltas)):

    val = p_vals[i]

    p_vals[i] -= deltas[i] / 2.
    mp.set_param_vals(p_vals)
    rev_uc = mp.get_model().get_unit_cell().parameters()
    rev_vec = mp.get_model().get_real_space_vectors()
    rev_vec = [Ut * vec for vec in rev_vec]
    rev_B = matrix.sqr(mp.get_model().get_B())
    rev_O = rev_B.transpose().inverse()

    p_vals[i] += deltas[i]
    mp.set_param_vals(p_vals)
    fwd_uc = mp.get_model().get_unit_cell().parameters()
    fwd_vec = mp.get_model().get_real_space_vectors()
    fwd_vec = [Ut * vec for vec in fwd_vec]
    fwd_B = matrix.sqr(mp.get_model().get_B())
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

from scitbx.math import angle_derivative_wrt_vectors

B = matrix.sqr(crystal.get_B())
O = (B.transpose()).inverse()
a, b, c, aa, bb, cc = crystal.get_unit_cell().parameters()
aa *= DEG2RAD
bb *= DEG2RAD
cc *= DEG2RAD
Ut = matrix.sqr(crystal.get_U()).transpose()
avec, bvec, cvec = [Ut * vec for vec in crystal.get_real_space_vectors()]

# calculate d[B^T]/dp
dB_dp = xluc_param.get_ds_dp()
dBT_dp = [dB.transpose() for dB in dB_dp]

# calculate d[O]/dp
dO_dp = [-O * dBT * O for dBT in dBT_dp]

# function to get analytical derivative of angles wrt vectors
def dangle(u, v):
  return [matrix.col(e) for e in angle_derivative_wrt_vectors(u,v)]

dalpha_db, dalpha_dc = dangle(bvec, cvec)
dbeta_da, dbeta_dc = dangle(avec, cvec)
dgamma_da, dgamma_db = dangle(avec, bvec)

# get all FD derivatives
fd_grad = check_fd_gradients(xluc_param)

# look at each parameter
for i, dO in enumerate(dO_dp):

  #print
  #print "***** PARAMETER {0} *****".format(i)

  #print "dB_dp analytical"
  #print dB_dp[i]
  #print "dB_dp FD"
  #print fd_grad[i]['dB_dp']
  #print

  # dB_dp is good. What about dO_dp?

  #print "O MATRIX"
  #print "dO_dp analytical"
  #print dO.round(6)
  #print "dO_dp FD"
  #print fd_grad[i]['dO_dp'].round(6)
  #print
  assert approx_equal(dO, fd_grad[i]['dO_dp'])

  # extract derivatives of each unit cell vector wrt p
  dav_dp, dbv_dp, dcv_dp = dO.transpose().as_list_of_lists()
  dav_dp = matrix.col(dav_dp)
  dbv_dp = matrix.col(dbv_dp)
  dcv_dp = matrix.col(dcv_dp)

  # check these are correct vs FD
  #print "CELL VECTORS"
  #diff = dav_dp - fd_grad[i]['davec_dp']
  #print 2 * diff.length() / (dav_dp.length() + fd_grad[i]['davec_dp'].length()) * 100
  #print 'davec_dp analytical: {0:.5f} {1:.5f} {2:.5f}'.format(*dav_dp.elems)
  #print 'davec_dp finite diff: {0:.5f} {1:.5f} {2:.5f}'.format(*fd_grad[i]['davec_dp'].elems)
  assert approx_equal(dav_dp, fd_grad[i]['davec_dp'])

  #diff = dbv_dp - fd_grad[i]['dbvec_dp']
  #print 2 * diff.length() / (dbv_dp.length() + fd_grad[i]['dbvec_dp'].length()) * 100
  #print 'dbvec_dp analytical: {0:.5f} {1:.5f} {2:.5f}'.format(*dbv_dp.elems)
  #print 'dbvec_dp finite diff: {0:.5f} {1:.5f} {2:.5f}'.format(*fd_grad[i]['dbvec_dp'].elems)
  assert approx_equal(dbv_dp, fd_grad[i]['dbvec_dp'])

  #diff = dcv_dp - fd_grad[i]['dcvec_dp']
  #print 2 * diff.length() / (dcv_dp.length() + fd_grad[i]['dcvec_dp'].length()) * 100
  #print 'dcvec_dp analytical: {0:.5f} {1:.5f} {2:.5f}'.format(*dcv_dp.elems)
  #print 'dcvec_dp finite diff: {0:.5f} {1:.5f} {2:.5f}'.format(*fd_grad[i]['dcvec_dp'].elems)
  #print
  assert approx_equal(dcv_dp, fd_grad[i]['dcvec_dp'])

  #print "CELL LENGTHS"
  da_dp = 1./a * avec.dot(dav_dp)
  #print "d[a]/dp{2} analytical: {0:.5f} FD: {1:.5f}".format(da_dp, fd_grad[i]['da_dp'], i)
  assert approx_equal(da_dp, fd_grad[i]['da_dp'])

  db_dp = 1./b * bvec.dot(dbv_dp)
  #print "d[b]/dp{2} analytical: {0:.5f} FD: {1:.5f}".format(db_dp, fd_grad[i]['db_dp'], i)
  assert approx_equal(db_dp, fd_grad[i]['db_dp'])

  dc_dp = 1./c * cvec.dot(dcv_dp)
  #print "d[c]/dp{2} analytical: {0:.5f} FD: {1:.5f}".format(dc_dp, fd_grad[i]['dc_dp'], i)
  assert approx_equal(dc_dp, fd_grad[i]['dc_dp'])

  #print
  #print "CELL ANGLES"

  daa_dp = RAD2DEG * (dbv_dp.dot(dalpha_db) + dcv_dp.dot(dalpha_dc))
  dbb_dp = RAD2DEG * (dav_dp.dot(dbeta_da) + dcv_dp.dot(dbeta_dc))
  dcc_dp = RAD2DEG * (dav_dp.dot(dgamma_da) + dbv_dp.dot(dgamma_db))

  #print "d[alpha]/dp{2} analytical: {0:.5f} FD: {1:.5f}".format(daa_dp, fd_grad[i]['daa_dp'], i)
  #print "d[beta]/dp{2} analytical: {0:.5f} FD: {1:.5f}".format(dbb_dp, fd_grad[i]['dbb_dp'], i)
  #print "d[gamma]/dp{2} analytical: {0:.5f} FD: {1:.5f}".format(dcc_dp, fd_grad[i]['dcc_dp'], i)
  assert approx_equal(daa_dp, fd_grad[i]['daa_dp'])
  assert approx_equal(dbb_dp, fd_grad[i]['dbb_dp'])
  assert approx_equal(dcc_dp, fd_grad[i]['dcc_dp'])

print "OK"

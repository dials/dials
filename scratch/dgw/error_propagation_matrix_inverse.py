#!/usr/bin/env python
#
#
#  Copyright (C) (2015) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman.
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division
import math
from scitbx import matrix
from scitbx.array_family import flex
import random

"""Implementation of the propagation of errors formula for matrix inversion
given in Lefebvre et al. (1999) http://arxiv.org/abs/hep-ex/9909031. As in
that paper the analytical formula is tested against Monte Carlo simulation."""

def random_vector_close_to(vector):
  vec = matrix.col(vector)
  vec2 = vec.rotate_around_origin(matrix.col(
              (random.random(),
               random.random(),
               random.random())).normalize(),
               random.gauss(0, 5), deg = True)
  length_multiplier = max(0.5, random.gauss(1, 0.1))
  return vec2 * length_multiplier

def create_Bmat():
  """Create a random crystal model, in P1 for maximum generality. Return the
  B matrix of this crystal"""

  from dxtbx.model.crystal import crystal_model

  vecs = map(random_vector_close_to,
             [(20, 0, 0),
              (0, 20, 0),
              (0, 0, 20)])

  return crystal_model(*vecs, space_group_symbol = "P 1").get_B()

def create_sig_mat(mat, fraction_sd=0.01):
  """Create a matrix of errors of the elements of mat, where each error is a
  normal deviate with a sigma some fraction of the absolute element value"""

  vals = [random.gauss(0.0, abs(fraction_sd*e)) for e in mat]
  return matrix.sqr(vals)

def perturb_mat(mat, fraction_sd=0.01):
  """Perturb the elements of mat by normal deviates with sigmas of some
  fraction of the absolute element values."""

  return mat + create_sig_mat(mat, fraction_sd)

def cov(a, b):
  """Return the sample covariance of vectors a and b"""
  a = flex.double(a)
  b = flex.double(b)
  n = len(a)
  assert n == len(b)
  resid_a = a - flex.mean(a)
  resid_b = b - flex.mean(b)
  return flex.sum(resid_a*resid_b) / (n - 1)

def calc_monte_carlo_covariances(mats):
  """Given the sequence of matrices mats, calculate the variance-covariance
  matrix of the elements"""

  # check input
  assert all([e.is_square() for e in mats])
  n = mats[0].n_rows()
  assert all([e.n_rows() == n for e in mats])

  # create an empty var-cov matrix
  covmat = flex.double(flex.grid(n**2,n**2), 0.0)

  for i in range(covmat.all()[0]):
    for j in range(covmat.all()[1]):
      a = [m[i] for m in mats]
      b = [m[j] for m in mats]
      covmat[i,j] = cov(a,b)

  return covmat

def calc_monte_carlo_population_covariances(mats, mean_matrix):
  """Given the sequence of matrices mats, calculate the variance-covariance
  matrix of the elements, using the known mean values for each elt in mat
  to avoid the approximation implied in taking the sample covariance"""

  # check input
  assert all([e.is_square() for e in mats])
  n = mats[0].n_rows()
  assert all([e.n_rows() == n for e in mats])

  # create an empty var-cov matrix
  covmat = flex.double(flex.grid(n**2,n**2), 0.0)

  for i in range(covmat.all()[0]):
    for j in range(covmat.all()[1]):
      resid_a = flex.double([m[i] - mean_matrix[i] for m in mats])
      resid_b = flex.double([m[j] - mean_matrix[j] for m in mats])
      covmat[i,j] = flex.mean(resid_a*resid_b)

  return covmat

def calc_analytical_covariances(mat, cov_mat):
  """Implement analytical formula of Lefebvre et al. (1999)
  http://arxiv.org/abs/hep-ex/9909031 to calculate the covariances of elements
  of mat^-1, given the covariances of mat itself. This is not the most efficient
  way to approach the formula, but suitable for initial tests"""

  # initialise covariance matrix
  assert mat.is_square()
  n = mat.n_rows()

  # use flex for nice 2D indexing
  inv_mat = flex.double(mat.inverse())
  inv_mat.reshape(flex.grid(n, n))

  inv_cov_mat = flex.double(flex.grid(n**2,n**2), 0.0)
  for alpha in range(n):
    for beta in range(n):
      for a in range(n):
        for b in range(n):
          # cov(m^-1[alpha, beta], m^-1[a, b])
          elt = 0.0
          for i in range(n):
            for j in range(n):
              for k in range(n):
                for l in range(n):
                  # index into cov_mat after flattening mat
                  x = i * n + j
                  y = k * n + l
                  elt += inv_mat[alpha, i] * inv_mat[j, beta] * \
                       inv_mat[a, k] * inv_mat[l, b] * \
                       cov_mat[x, y]
          # index into inv_cov_mat after flattening inv_mat
          x = alpha * n + beta
          y = a * n + b
          inv_cov_mat[x, y] = elt
  return inv_cov_mat


def test_lefebvre():
  """Run the test presented in part 4 of the paper Lefebvre et al. (1999),
  http://arxiv.org/abs/hep-ex/9909031."""

  # Simple test following the paper. First MC sim
  mat = matrix.sqr((0.7,0.2,0.4,0.6))
  sig_mat = mat*0.01
  p_mats = [perturb_mat(mat) for i in range(10000)]
  inv_p_mats = [m.inverse() for m in p_mats]

  # Here use the more exact formula for population covariance as we know the
  # expected means
  cov_inv_mat_MC = calc_monte_carlo_population_covariances(inv_p_mats,
    mean_matrix=mat.inverse())

  # Now analytical formula
  cov_mat = flex.double(flex.grid(4, 4), 0.0)
  cov_mat.matrix_diagonal_set_in_place(flex.double([e**2 for e in sig_mat]))
  cov_inv_mat = calc_analytical_covariances(mat, cov_mat)

  frac = (cov_inv_mat_MC - cov_inv_mat) / cov_inv_mat_MC

  # assert all fractional errors are within 5%
  assert all([e < 0.05 for e in frac])

  return

def test_B_matrix():
  """Test errors in an inverted B matrix, when the errors in B are known (and
  are independent)"""

  Bmat = create_Bmat()
  inv_Bmat = Bmat.inverse()

  # Note that elements in the strict upper triangle of B are zero
  # so their perturbed values will be zero too and the errors in these elements
  # will also be zero (this is what we want)

  # calculate the covariance of elements of inv_B by Monte Carlo simulation,
  # assuming that each element of B has an independent normal error given by a
  # sigma of 1% of the element value
  perturbed_Bmats = [perturb_mat(Bmat, fraction_sd=0.01) for i in range(10000)]
  invBs = [m.inverse() for m in perturbed_Bmats]
  #cov_invB = calc_monte_carlo_covariances(invBs)
  cov_invB_MC = calc_monte_carlo_population_covariances(invBs, inv_Bmat)

  print "Monte Carlo"
  print cov_invB_MC.as_scitbx_matrix()

  # Now calculate using the analytical formula. First need the covariance
  # matrix of B itself. This is the diagonal matrix of errors applied in the
  # simulation.
  n = Bmat.n_rows()
  cov_B = flex.double(flex.grid(n**2, n**2), 0.0)
  sig_B = Bmat * 0.01
  cov_B.matrix_diagonal_set_in_place(flex.double([e**2 for e in sig_B]))

  # Now can use the analytical formula
  cov_invB = calc_analytical_covariances(Bmat, cov_B)

  print "Analytical"
  print cov_invB.as_scitbx_matrix()

if __name__ == '__main__':

  test_lefebvre()
  print "OK"

  test_B_matrix()
  print "OK"

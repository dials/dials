"""Test derivatives of the Hamilton, Rollett and Sparks least-squares target
(http://doi.org/10.1107/S0365110X65000233)"""

from __future__ import absolute_import, division, print_function

from libtbx.test_utils import approx_equal
from scitbx.random import poisson_distribution, variate

def av_I(I, w, g):
  '''Given observations I for some reflection, their weights w and their
  inverse scales g, calculate the best average intensity according to the
  formula in HRS1965'''

  nobs = len(I)
  assert len(w) == nobs
  assert len(g) == nobs

  numerator = sum([a*b*c for (a,b,c) in zip(w, g, I)])
  denominator = sum([a*b*b for (a,b) in zip(w, g)])

  assert denominator > 0
  return numerator / denominator

def dg_dp(g, iparam):
  '''Calculate dg/dp for a parameter of the model, specified by its index. In
  this trivial case we make the g values the parameters themselves. Therefore
  the dg_dp vector is a delta function'''

  result = [0] * len(g)
  assert iparam < len(g)
  result[iparam] = 1

  return result

def grad_av_I(I, w, g, iparam):
  '''Calculate the gradient of the best average intensity with respect to a
  parameter of the model, specified by its index.'''

  dg = dg_dp(g, iparam)

  u = sum([a*b*c for (a,b,c) in zip(w, g, I)])
  v = sum([a*b*b for (a,b) in zip(w, g)])

  du_dp = sum([a*b*c for (a,b,c) in zip(w, dg, I)])
  dv_dp = sum([2.*a*b*c for (a,b,c) in zip(w, g, dg)])

  result = (du_dp * v - dv_dp * u) / (v**2)
  return result

def fd_grad_av_I(I, w, g, iparam):
  '''Calculate the finite difference approximation to the gradient of the best
  average intensity with respect to a parameter of the model, specified by its
  index.'''

  p = list(g)
  delta = 1.e-7
  assert iparam < len(p)

  p[iparam] -= 0.5 * delta
  rev_state = av_I(I, w, p)

  p[iparam] += delta
  fwd_state = av_I(I, w, p)

  fd_grad = (fwd_state - rev_state) / delta

  return fd_grad

def residual(I, w, g, iparam):
  '''Calculate a residual of the HRS target'''

  wl = w[iparam]
  gl = g[iparam]
  Il = I[iparam]
  mrgI = av_I(I, w, g)

  return Il - gl * mrgI

def grad_r(I, w, g, iparam):
  '''Calculate the first derivative of the residual of the HRS target with
  respect to a parameter of the model, specified by its index.'''

  dg = dg_dp(g, iparam)
  dgl = dg[iparam]
  wl = w[iparam]
  gl = g[iparam]
  mrgI = av_I(I, w, g)

  term1 = -1. * mrgI * dgl
  term2a = sum([2.*a*b*c for (a,b,c) in zip(w, g, dg)]) * mrgI
  term2b = sum([a*b*c for (a,b,c) in zip(w, dg, I)])
  term2c = sum([a*b*b for (a,b) in zip(w, g)])

  result = term1 + gl * (term2a - term2b) / term2c
  return result

def fd_grad_r(I, w, g, iparam):
  '''Calculate the finite difference approximation to the first derivative of
  the residual of the HRS target with respect to a parameter of the model,
  specified by its index.'''

  p = list(g)
  delta = 1.e-7
  assert iparam < len(p)

  p[iparam] -= 0.5 * delta
  rev_state = residual(I, w, p, iparam)

  p[iparam] += delta
  fwd_state = residual(I, w, p, iparam)

  fd_grad = (fwd_state - rev_state) / delta

  return fd_grad

def test():
  # Test the expression in dials_regression/doc/notes/scaling/scaling.tex
  # (as of revision 1537) for d<Ih>/dp. Simulate 10 measurements of a
  # reflection with different scales. Here the scale increases linearly between
  # each equally-spaced measurement, but the form of the scale variation
  # doesn't affect the derivative calculation.
  nobs = 10

  # known scale factors
  K = [1 + e/40. for e in range(nobs)]
  g = [1./e for e in K]

  # intensities
  means = [100 * e for e in K]
  I = [variate(poisson_distribution(e))() for e in means]

  # weights (inverse variances of I)
  w = [1./e for e in means]

  mrgI = av_I(I, w, g)

  for iparam in range(nobs):
    dmrgI_dp = grad_av_I(I, w, g, iparam)
    fd_dmrgI_dp = fd_grad_av_I(I, w, g, iparam)

    assert approx_equal(dmrgI_dp, fd_dmrgI_dp)

  # Now test the complete expression for the first derivative of the residuals
  # of the HRS target.
  for iparam in range(nobs):
    dr_dp = grad_r(I, w, g, iparam)
    fd_dr_dp = fd_grad_r(I, w, g, iparam)

    assert approx_equal(dr_dp, fd_dr_dp)

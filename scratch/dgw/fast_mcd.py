#!/usr/bin/env python

"""Testing functions for multivariate outlier rejection by the FAST-MCD
algorithm"""

# Want implementation of Mahalanobis distance that matches this R session:

#> x1 <- round(rnorm(10,3), 3)
#> x2 <- round(x1 + rnorm(10), 3)
#> x3 <- round(x2 + runif(10), 3)
#> x1
# [1] 3.853 2.401 2.253 3.067 1.887 3.293 3.995 2.559 2.785 2.228
#> x2
# [1] 4.294 1.915 1.315 4.641 1.611 2.838 3.696 1.337 2.853 2.434
#> x3
# [1] 4.785 2.352 2.023 4.978 2.329 3.101 4.494 2.204 3.468 3.075
#> obs <- cbind(x1, x2, x3)
#> S <- var(obs)
#> S
#          x1        x2       x3
#x1 0.5020374 0.6667232 0.633355
#x2 0.6667232 1.4434718 1.326026
#x3 0.6333550 1.3260262 1.248315
#> mahalanobis(obs, c(mean(x1), mean(x2), mean(x3)), S)
# [1] 2.1838336 1.9673401 1.3335029 4.9191627 2.1246818 5.3297995 4.9022487
# [8] 2.5335913 0.1952562 1.5105832

from scitbx.array_family import flex

def sample_covariance(a, b):
  """Calculate sample covariance of two vectors"""

  N = len(a)
  assert len(b) == N
  return flex.sum((a - flex.mean(a)) * (b - flex.mean(b))) / (N - 1)

x1 = flex.double((3.853, 2.401, 2.253, 3.067, 1.887, 3.293, 3.995, 2.559, 2.785, 2.228))
x2 = flex.double((4.294, 1.915, 1.315, 4.641, 1.611, 2.838, 3.696, 1.337, 2.853, 2.434))
x3 = flex.double((4.785, 2.352, 2.023, 4.978, 2.329, 3.101, 4.494, 2.204, 3.468, 3.075))


def cov(*args):
  """Calculate covariance matrix between the arguments (should be flex.double
  arrays of equal length)"""

  lens = [len(e) for e in args]
  assert all([e == lens[0] for e in lens])

  ncols = len(args)
  cov = flex.double(flex.grid(ncols, ncols))
  for i in range(ncols):
    for j in range(i, ncols):
      cov[i,j] = sample_covariance(args[i], args[j])

  cov.matrix_copy_upper_to_lower_triangle_in_place()
  return cov

# observation matrix (shifted by means) and its transpose
obs = flex.double((flex.grid(len(x1), 3)))
obs.matrix_paste_column_in_place(x1 - flex.mean(x1), 0)
obs.matrix_paste_column_in_place(x2 - flex.mean(x2), 1)
obs.matrix_paste_column_in_place(x3 - flex.mean(x3), 2)
obs_t = obs.matrix_transpose()

# covariance matrix
covmat = cov(x1, x2, x3)

# test Mahalanobis distance. This is an inefficient approach to calculating this
# because we do a full matrix product when we only want the diagonal elements of
# the result. More efficient approaches exist. See e.g.
# http://blogs.sas.com/content/iml/2012/02/22/how-to-compute-mahalanobis-distance-in-sas.html
maha=((obs.matrix_multiply(covmat.matrix_inversion())).matrix_multiply(obs_t)).matrix_diagonal()

# NB This is actually the squared Mahalanobis distance!

from libtbx.test_utils import approx_equal
R_result = [2.1838336, 1.9673401, 1.3335029, 4.9191627, 2.1246818,
            5.3297995, 4.9022487, 2.5335913, 0.1952562, 1.5105832]
assert approx_equal(list(maha), R_result)

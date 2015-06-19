#!/usr/bin/env python

"""Testing functions for multivariate outlier rejection by the FAST-MCD
algorithm"""

# Want implementation of Mahalanobis distance that matches the test in R:

#> ma <- cbind(1:6, 1:3)
#> ma
#     [,1] [,2]
#[1,]    1    1
#[2,]    2    2
#[3,]    3    3
#[4,]    4    1
#[5,]    5    2
#[6,]    6    3
#> (S <-  var(ma))
#     [,1] [,2]
#[1,]  3.5  0.8
#[2,]  0.8  0.8
#> ?var
#> mahalanobis(c(0, 0), 1:2, S)
#[1] 5.37037

from scitbx.array_family import flex

def sample_covariance(a, b):
  """Calculate sample covariance of two vectors"""

  N = len(a)
  assert len(b) == N
  return flex.sum((a - flex.mean(a)) * (b - flex.mean(b))) / (N - 1)

a = flex.double((1, 2, 3, 4, 5, 6))
b = flex.double((1, 2, 3, 1, 2, 3))

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

# observation matrix and its transpose
# FIXME probably better to create the output matrix once and use paste column in
# place to avoid excessive copying when ndim > 2
obs_t = flex.double.concatenate(a,b)
obs_t.reshape(flex.grid(2,6))
obs = obs_t.matrix_transpose()

# covariance matrix
covmat = cov(a,b)

# the R test looks at a single observation of (0,0) with centres of (1,2)
vec = flex.double((0-1, 0-2))
vec.reshape(flex.grid(1,2))
vec_t = vec.matrix_transpose()

# test Mahalanobis distance
maha=(vec.matrix_multiply(covmat.matrix_inversion())).matrix_multiply(vec_t)
assert round(maha[0],6) == 5.37037



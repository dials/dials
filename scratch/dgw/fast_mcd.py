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
from math import floor

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

def maha_dist_sq(cols, center, cov):
  """Calculate squared Mahalanobis distance of all observations (rows in the
  vectors contained in the list cols) from the center vector with respect to
  the covariance matrix cov"""

  n = len(cols[0])
  p = len(cols)
  assert len(center) == p

  # observation matrix (shifted by center) and its transpose
  obs = flex.double(flex.grid(n, p))
  for i, col in enumerate(cols):
    obs.matrix_paste_column_in_place(col - center[i], i)
  obs_t = obs.matrix_transpose()

  d2 = ((obs.matrix_multiply(cov.matrix_inversion())).matrix_multiply(obs_t)).matrix_diagonal()
  return d2

# test Mahalanobis distance. This is an inefficient approach to calculating this
# because we do a full matrix product when we only want the diagonal elements of
# the result. More efficient approaches exist. See e.g.
# http://blogs.sas.com/content/iml/2012/02/22/how-to-compute-mahalanobis-distance-in-sas.html
cols = [x1, x2, x3]
center = [flex.mean(e) for e in cols]
covmat = cov(x1, x2, x3)
maha = maha_dist_sq(cols, center, covmat)

# NB This is actually the squared Mahalanobis distance!

from libtbx.test_utils import approx_equal
R_result = [2.1838336, 1.9673401, 1.3335029, 4.9191627, 2.1246818,
            5.3297995, 4.9022487, 2.5335913, 0.1952562, 1.5105832]
assert approx_equal(list(maha), R_result)

# Now trial implementation of FAST-MCD algorithm

class FastMCD(object):
  """Experimental implementation of the FAST-MCD algorithm of Rousseeuw and
  van Driessen"""

  def __init__(self, data, alpha=0.5):
    """data expected to be a list of flex.double arrays of the same length,
    representing the vectors of observations in each dimension"""

    # the full dataset as separate vectors
    self._data = data

    # number of variables
    self._p = len(self._data)
    # p == 1 is the univariate case, best dealt with using a different (exact)
    # algorithm. TODO: implement that here too.
    assert self._p > 1

    # number of observations
    lens = [len(e) for e in self._data]
    assert all([e == lens[0] for e in lens])
    self._n = lens[0]

    # some input checks
    assert self._n > self._p

    # default initial subset size

    n2 = (self._n + self._p + 1) // 2
    self._h = int(floor(2 * n2 - self._n + 2 * (self._n - n2) *  alpha))
    # In the original FAST-MCD, if h == n it reports a single
    # location and scatter estimate for the whole dataset and stops. Currently
    # limit this implementation to h < n
    assert self._h < self._n

    # hardcoded for now...
    self._max_n_groups = 5
    self._min_group_size = 300
    self._n_trials = 500

    # iteration limits
    self._k1 = 2
    self._k2 = 2
    self._k3 = 100

  @staticmethod
  def means_and_covariance(vecs):
    """Prepare a dataset of equal length vectors for Mahalanobis distance
    squared calculation. The maha_dist_sq function requires the vectors,
    the vector of their means and their covariance matrix. Given the vectors,
    return the latter pair as a tuple"""

    center = [flex.mean(e) for e in vecs]
    covmat = cov(*vecs)
    return (center, covmat)

  @staticmethod
  def sample_data(data, sample_size):
    """sample (without replacement) the data vectors to select the same
    sample_size rows from each."""

    n = len(data[0])
    rows = flex.random_selection(n, sample_size)
    cols = [e.select(rows) for e in data]
    return cols

  def sample_data_and_group(self, data, sample_size, ngroups):
    """as sample_data, but follow by splitting each sampled vector into
    groups of approximately equal size."""

    sampled = self.sample_data(sample_size = sample_size)

    # random permutation
    p = flex.random_permutation(sample_size)
    permuted = [col.select(p) for col in data]

    # determine groups
    blocksize = int(sample_size / ngroups)
    rem = sample_size % ngroups
    blocksizes = [blocksize] * (ngroups - rem) + [blocksize + 1] * rem

    starts = [0]
    ends = [blocksizes[0]]
    for b in blocksizes[1:]:
      starts.append(ends[-1])
      ends.append(ends[-1] + b)
    blocks = zip(starts, ends)

    # split into groups
    groups = []
    for s,e in blocks:
      groups.append([col[s:e] for col in permuted])

    return groups

  def form_initial_subset(self, h, data):
    """Method 2 of subsection 3.1 of R&vD"""

    # draw random p+1 subset J (or larger if required)
    detS0 = 0.0
    i = 0
    while not detS0 > 0.0:
      J = self.sample_data(data, sample_size=self._p+1+i)
      i += 1
      T0, S0 = self.means_and_covariance(J)
      detS0 = S0.matrix_determinant_via_lu()

    H1 = self.concentration_step(h, data, T0, S0)
    return H1

  @staticmethod
  def concentration_step(h, data, T, S):
    """Practical application of Theorem 1 of R&vD"""

    d2s = maha_dist_sq(data, T, S)
    p = flex.sort_permutation(d2s)
    H1 = [col.select(p)[0:h] for col in data]
    return H1

  def small_dataset_estimate(self):
    """When a dataset is small, perform the initial trials directly on the
    whole dataset"""

    trials = []
    for i in xrange(self._n_trials):

      H1 = self.form_initial_subset(h=self._h, data=self._data)
      T1, S1 = self.means_and_covariance(H1)
      detS1 = S1.matrix_determinant_via_lu()

      # perform concentration steps
      detScurr, Tcurr, Scurr = detS1, T1, S1
      for j in xrange(self._k1): # take maximum of k1 steps

        Hnew = self.concentration_step(self._h, self._data, Tcurr, Scurr)
        Tnew, Snew = self.means_and_covariance(Hnew)
        detSnew = Snew.matrix_determinant_via_lu()

        # detS3 < detS2 < detS1 by Theorem 1. In practice (rounding errors?)
        # this is not always the case here. Ensure that detScurr is no smaller than
        # one billionth the value of detSnew less than detSnew
        assert detScurr > (detSnew - detSnew/1.e9)
        detScurr, Tcurr, Scurr = detSnew, Tnew, Snew

      trials.append((detSnew, Tnew, Snew))

    # choose 10 trials with the lowest detS3
    trials.sort(key=lambda x: x[0])
    best_trials = []
    for i in xrange(10):
      detCurr, Tcurr, Scurr = trials[i]
      for j in xrange(self._k3): # take maximum of k3 steps
        Hnew = self.concentration_step(self._h, self._data, Tcurr, Scurr)
        Tnew, Snew = self.means_and_covariance(Hnew)
        detNew = Snew.matrix_determinant_via_lu()
        if detNew == detCurr:
          print "trial {0}; iteration {1}; convergence".format(i,j)
          break
        detCurr, Tcurr, Scurr = detNew, Tnew, Snew
      best_trials.append((detCurr, Tnew, Snew))

    # Find the minimum covariance determinant from that set of 10
    best_trials.sort(key=lambda x: x[0])
    _, Tbest, Sbest = best_trials[0]
    return Tbest, Sbest

  def large_dataset_estimate(self):
    """When a dataset is large, construct disjoint subsets of the full data
    and perform initial trials within each of these then merge"""

    ngroups = int(self._n / self._min_group_size)
    if ngroups < self._max_n_groups:
      # use all the data and split into ngroups
      groups = sample_data_and_group(data=self._data,
                                     sample_size=self._n,
                                     ngroups=ngroups)
    else:
      # sample the data and split into the maximum number of groups
      sample_size = self._min_group_size * self._max_n_groups
      groups = sample_data_and_group(data=self._data,
                                     sample_size=sample_size,
                                     ngroups=self._max_n_groups)

    #### FIXME #####
    return Tbest, Sbest

  def detect_outliers(self):
    """Use the FAST-MCD estimate of location and scatter, T and S, to classify
    the rows of the input data as inlier or outlier. Report as a flex.bool"""

    # algorithm for a small number of observations (up to twice the minimum
    # group size)
    if self._n < 2 * self._min_group_size:
      T, S = self.small_dataset_estimate()

    # algorithm for a larger number of observations
    else:
      T, S = self.large_dataset_estimate()

    print T
    print S.as_scitbx_matrix()
    #vecs, center, covmat = self.means_and_covariance(self._data)
    #dists = maha_sq(vecs, T, S)

    #TODO use dists to classify as outliers
    return


# some test data, from R package robustbase: Hawkins, Bradu, Kass's Artificial Data
hbk = """10.1 19.6 28.3
 9.5 20.5 28.9
10.7 20.2 31.0
 9.9 21.5 31.7
10.3 21.1 31.1
10.8 20.4 29.2
10.5 20.9 29.1
 9.9 19.6 28.8
 9.7 20.7 31.0
 9.3 19.7 30.3
11.0 24.0 35.0
12.0 23.0 37.0
12.0 26.0 34.0
11.0 34.0 34.0
 3.4  2.9  2.1
 3.1  2.2  0.3
 0.0  1.6  0.2
 2.3  1.6  2.0
 0.8  2.9  1.6
 3.1  3.4  2.2
 2.6  2.2  1.9
 0.4  3.2  1.9
 2.0  2.3  0.8
 1.3  2.3  0.5
 1.0  0.0  0.4
 0.9  3.3  2.5
 3.3  2.5  2.9
 1.8  0.8  2.0
 1.2  0.9  0.8
 1.2  0.7  3.4
 3.1  1.4  1.0
 0.5  2.4  0.3
 1.5  3.1  1.5
 0.4  0.0  0.7
 3.1  2.4  3.0
 1.1  2.2  2.7
 0.1  3.0  2.6
 1.5  1.2  0.2
 2.1  0.0  1.2
 0.5  2.0  1.2
 3.4  1.6  2.9
 0.3  1.0  2.7
 0.1  3.3  0.9
 1.8  0.5  3.2
 1.9  0.1  0.6
 1.8  0.5  3.0
 3.0  0.1  0.8
 3.1  1.6  3.0
 3.1  2.5  1.9
 2.1  2.8  2.9
 2.3  1.5  0.4
 3.3  0.6  1.2
 0.3  0.4  3.3
 1.1  3.0  0.3
 0.5  2.4  0.9
 1.8  3.2  0.9
 1.8  0.7  0.7
 2.4  3.4  1.5
 1.6  2.1  3.0
 0.3  1.5  3.3
 0.4  3.4  3.0
 0.9  0.1  0.3
 1.1  2.7  0.2
 2.8  3.0  2.9
 2.0  0.7  2.7
 0.2  1.8  0.8
 1.6  2.0  1.2
 0.1  0.0  1.1
 2.0  0.6  0.3
 1.0  2.2  2.9
 2.2  2.5  2.3
 0.6  2.0  1.5
 0.3  1.7  2.2
 0.0  2.2  1.6
 0.3  0.4  2.6"""

# unpack the data into vectors
rows = [[float(e) for e in row.split()] for row in hbk.splitlines()]
x1, x2, x3 = [flex.double(e) for e in zip(*rows)]

fast_mcd = FastMCD([x1, x2, x3])
fast_mcd.detect_outliers()

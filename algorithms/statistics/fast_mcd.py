from __future__ import absolute_import, division, print_function

import math

from dials_refinement_helpers_ext import maha_dist_sq as maha_dist_sq_cpp
from dials_refinement_helpers_ext import mcd_consistency
from scitbx.array_family import flex


def sample_covariance(a, b):
    """Calculate sample covariance of two vectors"""

    N = len(a)
    assert len(b) == N
    return flex.sum((a - flex.mean(a)) * (b - flex.mean(b))) / (N - 1)


def cov(*args):
    """Calculate covariance matrix between the arguments (should be flex.double
    arrays of equal length)"""

    lens = [len(e) for e in args]
    assert all(e == lens[0] for e in lens)

    ncols = len(args)
    cov = flex.double(flex.grid(ncols, ncols))
    for i in range(ncols):
        for j in range(i, ncols):
            cov[i, j] = sample_covariance(args[i], args[j])

    cov.matrix_copy_upper_to_lower_triangle_in_place()
    return cov


def maha_dist_sq(cols, center, cov):
    """Calculate squared Mahalanobis distance of all observations (rows in the
    vectors contained in the list cols) from the center vector with respect to
    the covariance matrix cov"""

    n = len(cols[0])
    p = len(cols)
    assert len(center) == p

    # observation matrix
    obs = flex.double(flex.grid(n, p))
    for i, col in enumerate(cols):
        obs.matrix_paste_column_in_place(col, i)

    d2 = maha_dist_sq_cpp(obs, flex.double(center), cov)
    return d2


def mcd_finite_sample(p, n, alpha):
    """Finite sample correction factor for the MCD estimate. Described in
    Pison et al. Metrika (2002). doi.org/10.1007/s001840200191. Implementation
    based on 'rawcorfactor' in fastmcd.m from Continuous Sound and Vibration
    Analysis by Edward Zechmann"""

    from scitbx.lstbx import normal_eqns

    if p > 2:
        coeffqpkwad500 = [
            [-1.42764571687802, 1.26263336932151, 2],
            [-1.06141115981725, 1.28907991440387, 3],
        ]
        coeffqpkwad875 = [
            [-0.455179464070565, 1.11192541278794, 2],
            [-0.294241208320834, 1.09649329149811, 3],
        ]

        y_500 = [
            math.log(-coeffqpkwad500[0][0] / p ** coeffqpkwad500[0][1]),
            math.log(-coeffqpkwad500[1][0] / p ** coeffqpkwad500[1][1]),
        ]
        y_875 = [
            math.log(-coeffqpkwad875[0][0] / p ** coeffqpkwad875[0][1]),
            math.log(-coeffqpkwad875[1][0] / p ** coeffqpkwad875[1][1]),
        ]
        A_500 = [
            [1, -math.log(coeffqpkwad500[0][2] * p ** 2)],
            [1, -math.log(coeffqpkwad500[1][2] * p ** 2)],
        ]
        A_875 = [
            [1, -math.log(coeffqpkwad875[0][2] * p ** 2)],
            [1, -math.log(coeffqpkwad875[1][2] * p ** 2)],
        ]

        # solve the set of equations labelled _500
        eqs = normal_eqns.linear_ls(2)
        for i in range(2):
            eqs.add_equation(
                right_hand_side=y_500[i],
                design_matrix_row=flex.double(A_500[i]),
                weight=1,
            )
        eqs.solve()
        coeffic_500 = eqs.solution()

        # solve the set of equations labelled _875
        eqs = normal_eqns.linear_ls(2)
        for i in range(2):
            eqs.add_equation(
                right_hand_side=y_875[i],
                design_matrix_row=flex.double(A_875[i]),
                weight=1,
            )
        eqs.solve()
        coeffic_875 = eqs.solution()

        fp_500_n = 1 - math.exp(coeffic_500[0]) / n ** coeffic_500[1]
        fp_875_n = 1 - math.exp(coeffic_875[0]) / n ** coeffic_875[1]

    elif p == 2:
        fp_500_n = 1 - math.exp(0.673292623522027) / n ** 0.691365864961895
        fp_875_n = 1 - math.exp(0.446537815635445) / n ** 1.06690782995919

    elif p == 1:
        fp_500_n = 1 - math.exp(0.262024211897096) / n ** 0.604756680630497
        fp_875_n = 1 - math.exp(-0.351584646688712) / n ** 1.01646567502486

    if alpha <= 0.875:
        fp_alpha_n = fp_500_n + (fp_875_n - fp_500_n) / 0.375 * (alpha - 0.5)

    if 0.875 < alpha and alpha <= 1:
        fp_alpha_n = fp_875_n + (1 - fp_875_n) / 0.125 * (alpha - 0.875)

    return 1 / fp_alpha_n


class FastMCD(object):
    """Experimental implementation of the FAST-MCD algorithm of Rousseeuw and
    van Driessen"""

    def __init__(
        self,
        data,
        alpha=0.5,
        max_n_groups=5,
        min_group_size=300,
        n_trials=500,
        k1=2,
        k2=2,
        k3=100,
    ):
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
        assert all(e == lens[0] for e in lens)
        self._n = lens[0]

        # some input checks
        assert self._n > self._p

        # default initial subset size
        self._alpha = alpha
        n2 = (self._n + self._p + 1) // 2
        self._h = int(math.floor(2 * n2 - self._n + 2 * (self._n - n2) * alpha))
        # In the original FAST-MCD, if h == n it reports a single
        # location and scatter estimate for the whole dataset and stops. Currently
        # limit this implementation to h < n
        assert self._h < self._n

        # set sizes of groups and how many trials to perform in each
        self._max_n_groups = max_n_groups
        self._min_group_size = min_group_size
        self._n_trials = n_trials

        # iteration limits
        self._k1 = k1
        self._k2 = k2
        self._k3 = k3

        # correction factors
        self._consistency_fac = mcd_consistency(self._p, self._h / self._n)
        self._finite_samp_fac = mcd_finite_sample(self._p, self._n, self._alpha)

        # perform calculation
        self._T_raw = None
        self._S_raw = None

        self.run()

    def run(self):
        """Run the Fast MCD calculation"""

        # algorithm for a small number of observations (up to twice the minimum
        # group size)
        if self._n < 2 * self._min_group_size:
            self._T_raw, self._S_raw = self.small_dataset_estimate()

        # algorithm for a larger number of observations
        else:
            self._T_raw, self._S_raw = self.large_dataset_estimate()

    def get_raw_T_and_S(self):
        """Get the raw MCD location (T) and covariance matrix (S) estimates"""

        return self._T_raw, self._S_raw

    def get_corrected_T_and_S(self):
        """Get the MCD location (T) and covariance matrix (S) estimates corrected
        for normal model consistency and finite-sample size"""

        fac = self._consistency_fac * self._finite_samp_fac
        return self._T_raw, self._S_raw * fac

    @staticmethod
    def means_and_covariance(vecs):
        """Prepare a dataset of equal length vectors for Mahalanobis distance
        squared calculation. The maha_dist_sq function requires the vectors,
        the vector of their means and their covariance matrix. Given the vectors,
        return the latter pair as a tuple"""

        center = flex.double([flex.mean(e) for e in vecs])
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

    def split_into_groups(self, sample, ngroups):
        """Split each vector in the data sample into groups of approximately equal
        size."""

        # number of obs in the sample
        sample_size = len(sample[0])

        # random permutation
        p = flex.random_permutation(sample_size)
        permuted = [col.select(p) for col in sample]

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
        groups = [[col[start:end] for col in permuted] for start, end in blocks]

        return groups

    def form_initial_subset(self, h, data):
        """Method 2 of subsection 3.1 of R&vD"""

        # permutation of input data for sampling
        p = flex.random_permutation(len(data[0]))
        permuted = [col.select(p) for col in data]

        # draw random p+1 subset J (or larger if required)
        detS0 = 0.0
        i = 0
        while not detS0 > 0.0:
            subset_size = self._p + 1 + i
            J = [e[0:subset_size] for e in permuted]
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
        for i in range(self._n_trials):

            H1 = self.form_initial_subset(h=self._h, data=self._data)
            T1, S1 = self.means_and_covariance(H1)
            detS1 = S1.matrix_determinant_via_lu()

            # perform concentration steps
            detScurr, Tcurr, Scurr = detS1, T1, S1
            for j in range(self._k1):  # take maximum of k1 steps

                Hnew = self.concentration_step(self._h, self._data, Tcurr, Scurr)
                Tnew, Snew = self.means_and_covariance(Hnew)
                detSnew = Snew.matrix_determinant_via_lu()

                # detS3 < detS2 < detS1 by Theorem 1. In practice (rounding errors?)
                # this is not always the case here. Ensure that detScurr is no smaller than
                # one billionth the value of detSnew less than detSnew
                assert detScurr > (detSnew - detSnew / 1.0e9)
                detScurr, Tcurr, Scurr = detSnew, Tnew, Snew

            trials.append((detSnew, Tnew, Snew))

        # choose 10 trials with the lowest detS3
        trials.sort(key=lambda x: x[0])
        best_trials = []
        for i in range(10):
            detCurr, Tcurr, Scurr = trials[i]
            for j in range(self._k3):  # take maximum of k3 steps
                Hnew = self.concentration_step(self._h, self._data, Tcurr, Scurr)
                Tnew, Snew = self.means_and_covariance(Hnew)
                detNew = Snew.matrix_determinant_via_lu()
                if detNew == detCurr:
                    # print "trial {0}; iteration {1}; convergence".format(i,j)
                    break
                detCurr, Tcurr, Scurr = detNew, Tnew, Snew
            best_trials.append((detNew, Tnew, Snew))

        # Find the minimum covariance determinant from that set of 10
        best_trials.sort(key=lambda x: x[0])
        _, Tbest, Sbest = best_trials[0]
        return Tbest, Sbest

    def large_dataset_estimate(self):
        """When a dataset is large, construct disjoint subsets of the full data
        and perform initial trials within each of these, then merge"""

        ngroups = int(self._n / self._min_group_size)
        if ngroups < self._max_n_groups:
            # use all the data and split into ngroups
            sample_size = self._n
        else:
            # sample the data and split into the maximum number of groups
            ngroups = self._max_n_groups
            sample_size = self._min_group_size * self._max_n_groups

        # sample the data and split into groups
        sampled = self.sample_data(self._data, sample_size=sample_size)
        groups = self.split_into_groups(sample=sampled, ngroups=ngroups)

        # work within the groups now
        n_trials = self._n_trials // ngroups
        trials = []
        h_frac = self._h / self._n
        for group in groups:

            h_sub = int(len(group[0]) * h_frac)
            gp_trials = []
            for i in range(n_trials):

                H1 = self.form_initial_subset(h=h_sub, data=group)
                T1, S1 = self.means_and_covariance(H1)
                detS1 = S1.matrix_determinant_via_lu()

                # perform concentration steps
                detScurr, Tcurr, Scurr = detS1, T1, S1
                for j in range(self._k1):  # take k1 steps

                    Hnew = self.concentration_step(h_sub, group, Tcurr, Scurr)
                    Tnew, Snew = self.means_and_covariance(Hnew)
                    detSnew = Snew.matrix_determinant_via_lu()

                    # detS3 < detS2 < detS1 by Theorem 1. In practice (rounding errors?)
                    # this is not always the case here. Ensure that detScurr is no smaller than
                    # one billionth the value of detSnew less than detSnew
                    assert detScurr > (detSnew - detSnew / 1.0e9)
                    detScurr, Tcurr, Scurr = detSnew, Tnew, Snew

                gp_trials.append((detSnew, Tnew, Snew))

            # choose 10 trials with the lowest determinant and put in the outer list
            gp_trials.sort(key=lambda x: x[0])
            trials.extend(gp_trials[0:10])

        # now have 10 best trials from each group. Work with the merged (==sampled)
        # set
        mrgd_trials = []
        h_mrgd = int(sample_size * h_frac)
        for trial in trials:

            detScurr, Tcurr, Scurr = trial
            for j in range(self._k2):  # take k2 steps

                Hnew = self.concentration_step(h_mrgd, sampled, Tcurr, Scurr)
                Tnew, Snew = self.means_and_covariance(Hnew)
                detSnew = Snew.matrix_determinant_via_lu()
                detScurr, Tcurr, Scurr = detSnew, Tnew, Snew

            mrgd_trials.append((detSnew, Tnew, Snew))

        # sort trials by the lowest detS3 and work with the whole dataset now
        mrgd_trials.sort(key=lambda x: x[0])

        # choose number of steps to iterate based on dataset size (ugly)
        size = self._n * self._p
        if size <= 100000:
            k4 = self._k3
        elif size <= 200000:
            k4 = 10
        elif size <= 300000:
            k4 = 9
        elif size <= 400000:
            k4 = 8
        elif size <= 500000:
            k4 = 7
        elif size <= 600000:
            k4 = 6
        elif size <= 700000:
            k4 = 5
        elif size <= 800000:
            k4 = 4
        elif size <= 900000:
            k4 = 3
        elif size <= 1000000:
            k4 = 2
        else:
            k4 = 1

        # choose number of trials to look at based on number of obs (ugly)
        n_reps = 1 if self._n > 5000 else 10

        best_trials = []
        for i in range(n_reps):
            detCurr, Tcurr, Scurr = mrgd_trials[i]
            for j in range(k4):  # take maximum of k4 steps
                Hnew = self.concentration_step(self._h, self._data, Tcurr, Scurr)
                Tnew, Snew = self.means_and_covariance(Hnew)
                detNew = Snew.matrix_determinant_via_lu()
                if detNew == detCurr:
                    # print "trial {0}; iteration {1}; convergence".format(i,j)
                    break
                detCurr, Tcurr, Scurr = detNew, Tnew, Snew
            # print "determinant", detNew
            # print "location:", list(Tnew)
            # print "covariance:"
            # print Snew.as_scitbx_matrix()
            best_trials.append((detNew, Tnew, Snew))

        # Find the minimum covariance determinant from that set of 10
        best_trials.sort(key=lambda x: x[0])
        _, Tbest, Sbest = best_trials[0]
        return Tbest, Sbest

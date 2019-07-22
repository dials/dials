"""Algorithms for determination of Laue group symmetry."""
from __future__ import division, absolute_import, print_function

import logging

logger = logging.getLogger(__name__)

import math

import scipy.stats

import libtbx
from libtbx import table_utils
from scitbx.array_family import flex
from scitbx.math import five_number_summary

from cctbx import crystal, sgtbx

from dials.algorithms.symmetry import symmetry_base


class determine_space_group(symmetry_base):
    """Determination of Laue group symmetry using algorithms similar to POINTLESS.

    See also:
      `Evans, P. (2006). Acta Cryst. D62, 72-82
      <https://doi.org/10.1107/S0907444905036693>`_ and
      `Evans, P. R. (2011). Acta Cryst. D67, 282-292
      <https://doi.org/10.1107/S090744491003982X>`_.

    """

    def __init__(
        self,
        intensities,
        normalisation="ml_aniso",
        lattice_symmetry_max_delta=2.0,
        d_min=libtbx.Auto,
        min_i_mean_over_sigma_mean=4,
        min_cc_half=0.6,
        relative_length_tolerance=None,
        absolute_angle_tolerance=None,
    ):
        """Intialise a determine_space_group object.

        Args:
          intensities (cctbx.miller.array): The intensities on which to perform
            symmetry anaylsis.
          normalisation (str): The normalisation method to use. Possible choices are
            'kernel', 'quasi', 'ml_iso' and 'ml_aniso'. Set to None to switch off
            normalisation altogether.
          lattice_symmetry_max_delta (float): The maximum value of delta for
            determining the lattice symmetry using the algorithm of Le Page (1982).
          d_min (float): Optional resolution cutoff to be applied to the input
            intensities. If set to :data:`libtbx.Auto` then d_min will be
            automatically determined according to the parameters
            ``min_i_mean_over_sigma_mean`` and ``min_cc_half``.
          min_i_mean_over_sigma_mean (float): minimum value of |I|/|sigma(i)| for
            automatic determination of resolution cutoff.
          min_cc_half (float): minimum value of CC1/2 for automatic determination of
            resolution cutoff.
          relative_length_tolerance (float): Relative length tolerance in checking
            consistency of input unit cells against the median unit cell.
          absolute_angle_tolerance (float): Absolute angle tolerance in checking
            consistency of input unit cells against the median unit cell.

        """
        super(determine_space_group, self).__init__(
            intensities,
            normalisation=normalisation,
            lattice_symmetry_max_delta=lattice_symmetry_max_delta,
            d_min=d_min,
            min_i_mean_over_sigma_mean=min_i_mean_over_sigma_mean,
            min_cc_half=min_cc_half,
            relative_length_tolerance=relative_length_tolerance,
            absolute_angle_tolerance=absolute_angle_tolerance,
        )

        self._estimate_cc_sig_fac()
        self._estimate_cc_true()
        self._score_symmetry_elements()
        self._score_laue_groups()

    def _estimate_cc_sig_fac(self):
        """Estimation of sigma(CC) as a function of sample size.

        Estimate the error in the correlation coefficient, sigma(CC) by using
        pairs of reflections at similar resolutions that are not related by
        potential symmetry. Using pairs of unrelated reflections at similar
        resolutions, calculate sigma(CC) == rms(CC) for groups of size N = 3..200.
        The constant CCsigFac is obtained from a linear fit of
        sigma(CC) to 1/N^(1/2), i.e.:
            sigma(CC) = CCsigFac/N^(1/2)

        """

        max_bins = 500
        reflections_per_bin = max(
            200, int(math.ceil(self.intensities.size() / max_bins))
        )
        binner = self.intensities.setup_binner_counting_sorted(
            reflections_per_bin=reflections_per_bin
        )

        a = flex.double()
        b = flex.double()
        ma_tmp = self.intensities.customized_copy(
            crystal_symmetry=crystal.symmetry(
                space_group=self.lattice_group,
                unit_cell=self.intensities.unit_cell(),
                assert_is_compatible_unit_cell=False,
            )
        ).map_to_asu()
        for i in range(binner.n_bins_all()):
            count = binner.counts()[i]
            if count == 0:
                continue
            bin_isel = binner.array_indices(i)
            p = flex.random_permutation(count)
            p = p[: 2 * (count // 2)]  # ensure even count
            ma_a = ma_tmp.select(bin_isel.select(p[: count // 2]))
            ma_b = ma_tmp.select(bin_isel.select(p[count // 2 :]))
            # only choose pairs of reflections that don't have the same indices
            # in the asu of the lattice group
            sel = ma_a.indices() != ma_b.indices()
            a.extend(ma_a.data().select(sel))
            b.extend(ma_b.data().select(sel))

        perm = flex.random_selection(a.size(), min(20000, a.size()))
        a = a.select(perm)
        b = b.select(perm)

        self.corr_unrelated = CorrelationCoefficientAccumulator(a, b)

        n_pairs = a.size()
        min_num_groups = 10  # minimum number of groups
        max_n_group = int(min(n_pairs / min_num_groups, 200))  # maximum number in group
        min_n_group = int(min(5, max_n_group))  # minimum number in group

        mean_ccs = flex.double()
        rms_ccs = flex.double()
        ns = flex.double()
        for n in range(min_n_group, max_n_group):
            ns.append(n)
            ccs = flex.double()
            for i in range(200):
                isel = flex.random_selection(a.size(), n)
                corr = CorrelationCoefficientAccumulator(a.select(isel), b.select(isel))
                ccs.append(corr.coefficient())

            mean_ccs.append(flex.mean(ccs))
            rms_ccs.append(flex.mean(flex.pow2(ccs)) ** 0.5)

        x = 1 / flex.pow(ns, 0.5)
        y = rms_ccs
        fit = flex.linear_regression(x, y)

        assert fit.is_well_defined()
        self.cc_sig_fac = fit.slope()

    def _estimate_cc_true(self):

        # A1.2. Estimation of E(CC; S).

        # (i)

        var_intensities = flex.mean_and_variance(
            self.intensities.data()
        ).unweighted_sample_variance()
        var_sigmas = flex.mean_and_variance(flex.pow2(self.intensities.sigmas())).mean()
        self.E_cc_true = var_intensities / (var_intensities + var_sigmas)

        # (ii)

        reindexed_intensities = self.intensities.change_basis(
            sgtbx.change_of_basis_op("-x,-y,-z")
        ).map_to_asu()
        x, y = self.intensities.common_sets(
            reindexed_intensities, assert_is_similar_symmetry=False
        )
        self.cc_identity = CorrelationCoefficientAccumulator(x.data(), y.data())

        min_sd = 0.05
        min_sample = 10
        sigma_1 = max(min_sd, self.cc_sig_fac / 200 ** 0.5)
        w1 = 0
        w2 = 0
        if sigma_1 > 0.0001:
            w1 = 1 / sigma_1 ** 2
        if self.cc_identity.n() > min_sample:
            sigma_2 = max(min_sd, self.cc_sig_fac / self.cc_identity.n() ** 0.5)
            w2 = 1 / sigma_2 ** 2

        assert (w1 + w2) > 0
        self.cc_true = (w1 * self.E_cc_true + w2 * self.cc_identity.coefficient()) / (
            w1 + w2
        )

        logger.debug("cc_true = w1 * E_cc_true + w2 * cc_identity)/(w1 + w2)")
        logger.debug("w1: %g", w1)
        logger.debug("w2: %g", w2)
        logger.debug("E_cc_true: %g", self.E_cc_true)
        logger.debug("cc_identity: %g", self.cc_identity.coefficient())
        logger.debug("cc_true: %g", self.cc_true)

    def _score_symmetry_elements(self):
        self.sym_op_scores = []
        for smx in self.lattice_group.smx():
            if smx.r().info().sense() < 0:
                continue
            self.sym_op_scores.append(
                ScoreSymmetryElement(
                    self.intensities, smx, self.cc_true, self.cc_sig_fac
                )
            )

    def _score_laue_groups(self):
        subgroup_scores = [
            ScoreSubGroup(subgrp, self.sym_op_scores)
            for subgrp in self.subgroups.result_groups
        ]
        total_likelihood = sum(score.likelihood for score in subgroup_scores)
        for score in subgroup_scores:
            score.likelihood /= total_likelihood
        self.subgroup_scores = sorted(
            subgroup_scores, key=lambda score: score.likelihood, reverse=True
        )

        # The 'confidence' scores are derived from the total probability of the best
        # solution p_best and that for the next best solution p_next:
        #   confidence = [p_best * (p_best - p_next)]^1/2.

        for i, score in enumerate(self.subgroup_scores[:-1]):
            next_score = self.subgroup_scores[i + 1]
            if score.likelihood > 0 and next_score.likelihood > 0:
                lgc = score.likelihood * (score.likelihood - next_score.likelihood)
                confidence = abs(lgc) ** 0.5
                if lgc < 0:
                    confidence = -confidence
                score.confidence = confidence

        self.best_solution = self.subgroup_scores[0]

    def __str__(self):
        """Return a string representation of the results.

        Returns:
          str:

        """
        output = []
        output.append("Input crystal symmetry:")
        output.append(str(self.input_intensities[0].space_group_info()))
        output.append(str(self.median_unit_cell))
        output.append("Change of basis op to minimum cell: %s" % self.cb_op_inp_min)
        output.append("Crystal symmetry in minimum cell:")
        output.append(str(self.intensities.space_group_info()))
        output.append(str(self.intensities.unit_cell()))
        output.append("Lattice point group: %s" % self.lattice_group.info())
        output.append(
            "Overall CC for %i unrelated pairs: %.3f"
            % (self.corr_unrelated.n(), self.corr_unrelated.coefficient())
        )
        output.append(
            "Estimated expectation value of true correlation coefficient E(CC) = %.3f"
            % self.E_cc_true
        )
        output.append("Estimated sd(CC) = %.3f / sqrt(N)" % self.cc_sig_fac)
        output.append(
            "Estimated E(CC) of true correlation coefficient from identity = %.3f"
            % self.cc_true
        )

        header = ("likelihood", "Z-CC", "CC", "N", "", "Operator")
        rows = [header]
        for score in self.sym_op_scores:
            if score.likelihood > 0.9:
                stars = "***"
            elif score.likelihood > 0.7:
                stars = "**"
            elif score.likelihood > 0.5:
                stars = "*"
            else:
                stars = ""
            rows.append(
                (
                    "%.3f" % score.likelihood,
                    "%.2f" % score.z_cc,
                    "%.2f" % score.cc.coefficient(),
                    "%i" % score.n_refs,
                    stars,
                    "%s" % score.sym_op.r().info(),
                )
            )
        output.append("Scoring individual symmetry elements")
        output.append(table_utils.format(rows, has_header=True, delim="  "))

        header = (
            "Patterson group",
            "",
            "Likelihood",
            "NetZcc",
            "Zcc+",
            "Zcc-",
            "CC",
            "CC-",
            "delta",
            "Reindex operator",
        )
        rows = [header]
        for score in self.subgroup_scores:
            if score.likelihood > 0.8:
                stars = "***"
            elif score.likelihood > 0.6:
                stars = "**"
            elif score.likelihood > 0.4:
                stars = "*"
            else:
                stars = ""
            rows.append(
                (
                    "%s" % score.subgroup["best_subsym"].space_group_info(),
                    stars,
                    "%.3f" % score.likelihood,
                    "% .2f" % score.z_cc_net,
                    "% .2f" % score.z_cc_for,
                    "% .2f" % score.z_cc_against,
                    "% .2f" % score.cc_for.coefficient(),
                    "% .2f" % score.cc_against.coefficient(),
                    "%.1f" % score.subgroup["max_angular_difference"],
                    "%s" % (score.subgroup["cb_op_inp_best"]),
                )
            )
        output.append("Scoring all possible sub-groups")
        output.append(table_utils.format(rows, has_header=True, delim="  "))

        output.append(
            "Best solution: %s"
            % self.best_solution.subgroup["best_subsym"].space_group_info()
        )
        output.append(
            "Unit cell: %s" % self.best_solution.subgroup["best_subsym"].unit_cell()
        )
        output.append(
            "Reindex operator: %s" % (self.best_solution.subgroup["cb_op_inp_best"])
        )
        output.append("Laue group probability: %.3f" % self.best_solution.likelihood)
        output.append("Laue group confidence: %.3f" % self.best_solution.confidence)
        return "\n".join(output)

    def as_dict(self):
        """Return a dictionary representation of the results.

        Returns:
          dict

        """
        d = {
            "input_symmetry": {
                "hall_symbol": self.input_intensities[0]
                .space_group()
                .type()
                .hall_symbol(),
                "unit_cell": self.median_unit_cell.parameters(),
            },
            "cb_op_inp_min": self.cb_op_inp_min.as_xyz(),
            "min_cell_symmetry": {
                "hall_symbol": self.intensities.space_group().type().hall_symbol(),
                "unit_cell": self.intensities.unit_cell().parameters(),
            },
            "lattice_point_group": self.lattice_group.type().hall_symbol(),
            "cc_unrelated_pairs": self.corr_unrelated.coefficient(),
            "n_unrelated_pairs": self.corr_unrelated.n(),
            "E_cc_true": self.E_cc_true,
            "cc_sig_fac": self.cc_sig_fac,
            "cc_true": self.cc_true,
        }

        d["sym_op_scores"] = [score.as_dict() for score in self.sym_op_scores]
        d["subgroup_scores"] = [score.as_dict() for score in self.subgroup_scores]
        return d

    def as_json(self, filename=None, indent=2):
        """Return a json representation of the results.

        Args:
          filename (str): Optional filename to export the json representation of
            the results.
          indent (int): The indent level for pretty-printing of the json. If ``None``
            is the most compact representation.

        Returns:
          str:

        """
        d = self.as_dict()
        import json

        json_str = json.dumps(d, indent=indent)
        if filename is not None:
            with open(filename, "wb") as f:
                f.write(json_str)
        return json.dumps(d, indent=indent)


class ScoreCorrelationCoefficient(object):
    def __init__(self, cc, sigma_cc, expected_cc, lower_bound=-1, upper_bound=1, k=2):
        self.cc = cc
        self.sigma_cc = sigma_cc
        self.expected_cc = expected_cc
        self._lower_bound = lower_bound
        self._upper_bound = upper_bound
        self._k = k
        self._compute_p_cc_given_s()
        self._compute_p_cc_given_not_s()

    def _compute_p_cc_given_s(self):
        self._p_cc_given_s = trunccauchy_pdf(
            self.cc,
            self._lower_bound,
            self._upper_bound,
            self.expected_cc,
            self.sigma_cc,
        )

    def _compute_p_cc_given_not_s(self):
        sump = scipy.integrate.quad(self._numerator, 0, 1)[0]
        sumw = scipy.integrate.quad(self._denominator, 0, 1)[0]
        self._p_cc_given_not_s = sump / sumw

    @property
    def p_cc_given_s(self):
        """Probability of observing this CC if the sym op is present, p(CC; S).

        Modelled by a Cauchy distribution centred on cc_true and width gamma = sigma_cc

        """
        return self._p_cc_given_s

    @property
    def p_cc_given_not_s(self):
        """Probability of observing this CC if the sym op is NOT present, p(CC; !S).

        """
        return self._p_cc_given_not_s

    @property
    def p_s_given_cc(self):
        """The likelihood of this symmetry element being present.

        p(S; CC) = p(CC; S) / (p(CC; S) + p(CC; !S))

        """
        return self._p_cc_given_s / (self._p_cc_given_s + self._p_cc_given_not_s)

    def _p_mu_power_pdf(self, m):
        return (1.0 - pow(m, self._k)) ** (1.0 / self._k)

    def _numerator(self, x):
        return trunccauchy_pdf(
            self.cc, self._lower_bound, self._upper_bound, loc=x, scale=self.sigma_cc
        ) * self._p_mu_power_pdf(x)

    def _denominator(self, x):
        return self._p_mu_power_pdf(x)


class ScoreSymmetryElement(object):
    """Analyse intensities for presence of a given symmetry operation.

    1) Calculate the correlation coefficient, CC, for the given sym op.

    2) Calculate the probability of observing this CC if the sym op is present,
       p(CC; S), modelled by a Cauchy distribution centred on cc_true and width
       gamma = sigma_cc.

    3) Calculate the probability of observing this CC if the sym op is
       NOT present, p(CC; !S).

    4) Calculate the likelihood of symmetry element being present,
       p(S; CC) = p(CC; S) / (p(CC; S) + p(CC; !S))

    See appendix A1 of `Evans, P. R. (2011). Acta Cryst. D67, 282-292.
    <https://doi.org/10.1107/S090744491003982X>`_

    """

    def __init__(self, intensities, sym_op, cc_true, cc_sig_fac):
        """Initialise a ScoreSymmetryElement object.

        Args:
          intensities (cctbx.miller.array): The intensities on which to perform
            symmetry anaylsis.
          sym_op (cctbx.sgtbx.rt_mx): The symmetry operation for analysis.
          cc_true (float): the expected value of CC if the symmetry element is present,
            E(CC; S)
          cc_sig_fac (float): Estimation of sigma(CC) as a function of sample size.

        """
        self.sym_op = sym_op
        assert self.sym_op.r().info().sense() >= 0
        self.cc = CorrelationCoefficientAccumulator()
        cb_op = sgtbx.change_of_basis_op(self.sym_op)
        cb_ops = [cb_op]
        if self.sym_op.r().order() > 2:
            # include inverse symmetry operation
            cb_ops.append(cb_op.inverse())
        for cb_op in cb_ops:
            if cb_op.is_identity_op():
                cb_op = sgtbx.change_of_basis_op("-x,-y,-z")
            reindexed_intensities = intensities.change_basis(cb_op).map_to_asu()
            x, y = intensities.common_sets(
                reindexed_intensities, assert_is_similar_symmetry=False
            )
            sel = sgtbx.space_group().expand_smx(self.sym_op).epsilon(x.indices()) == 1
            x = x.select(sel)
            y = y.select(sel)

            outliers = flex.bool(len(x.data()), False)
            iqr_multiplier = 20  # very generous tolerance
            for col in (x.data(), y.data()):
                min_x, q1_x, med_x, q3_x, max_x = five_number_summary(col)
                iqr_x = q3_x - q1_x
                cut_x = iqr_multiplier * iqr_x
                outliers.set_selected(col > q3_x + cut_x, True)
                outliers.set_selected(col < q1_x - cut_x, True)
            if outliers.count(True):
                logger.debug(
                    "Rejecting %s outlier value%s"
                    % (libtbx.utils.plural_s(outliers.count(True)))
                )
                x = x.select(~outliers)
                y = y.select(~outliers)

            self.cc += CorrelationCoefficientAccumulator(x.data(), y.data())

        self.n_refs = self.cc.n()
        if self.n_refs <= 0:
            self.likelihood = 0
            self.z_cc = 0
            return

        self.sigma_cc = max(0.1, cc_sig_fac / self.n_refs ** 0.5)
        self.z_cc = self.cc.coefficient() / self.sigma_cc
        score_cc = ScoreCorrelationCoefficient(
            self.cc.coefficient(), self.sigma_cc, cc_true
        )
        self.p_cc_given_s = score_cc.p_cc_given_s
        self.p_cc_given_not_s = score_cc.p_cc_given_not_s
        self.likelihood = score_cc.p_s_given_cc

    def __str__(self):
        """Return a string representation of the symmetry element scoring.

        Returns:
          str:

        """
        return "%.3f %.2f %.2f %i %s" % (
            self.likelihood,
            self.z_cc,
            self.cc.coefficient(),
            self.n_refs,
            self.sym_op.r().info(),
        )

    def as_dict(self):
        """Return a dictionary representation of the symmetry element scoring.

        The dictionary will contain the following keys:
          - likelihood: The likelihood of the symmetry element being present
          - z_cc: The Z-score for the correlation coefficent
          - cc: The correlation coefficient for the symmetry element
          - n_ref: The number of reflections contributing to the correlation
            coefficient
          - operator: The xyz representation of the symmetry element

        Returns:
          dict:

        """
        return {
            "likelihood": self.likelihood,
            "z_cc": self.z_cc,
            "cc": self.cc.coefficient(),
            "n_ref": self.n_refs,
            "operator": self.sym_op.as_xyz(),
        }


class ScoreSubGroup(object):
    """Score the probability of a given subgroup being the true subgroup.

    1) Calculates the combined correlation coefficients for symmetry operations
       present/absent from the subgroup.

    2) Calculates overall Zcc scores for symmetry elements present/absent from
       the subgroup.

    3) Calculates the overall likelihood for this subgroup.

    See appendix A2 of `Evans, P. R. (2011). Acta Cryst. D67, 282-292.
    <https://doi.org/10.1107/S090744491003982X>`_

    """

    def __init__(self, subgroup, sym_op_scores):
        """Initialise a ScoreSubGroup object.

        Args:
          subgroup (dict): A dictionary describing the subgroup as generated by
            :class:`cctbx.sgtbx.lattice_symmetry.metric_subgroups`.
          sym_op_scores (list): A list of :class:`ScoreSymmetryElement` objects for each
            symmetry element possibly in the lattice symmetry.

        """
        # Combined correlation coefficients for symmetry operations
        # present/absent from subgroup
        self.subgroup = subgroup
        patterson_group = subgroup["subsym"].space_group()
        self.cc_for = CorrelationCoefficientAccumulator()
        self.cc_against = CorrelationCoefficientAccumulator()
        for score in sym_op_scores:
            if score.sym_op in patterson_group:
                self.cc_for += score.cc
            else:
                self.cc_against += score.cc

        # Overall Zcc scores for symmetry elements present/absent from subgroup
        self.z_cc_for = 0
        self.z_cc_against = 0
        n_for = 0
        n_against = 0
        PL_for = 0
        PL_against = 0
        power = 2
        for score in sym_op_scores:
            if score.n_refs <= 2:
                continue
            if score.sym_op in patterson_group:
                self.z_cc_for += score.z_cc ** power
                n_for += 1
                PL_for += math.log(score.p_cc_given_s)
            else:
                self.z_cc_against += score.z_cc ** power
                n_against += 1
                PL_against += math.log(score.p_cc_given_not_s)

        # Overall likelihood for this subgroup
        self.likelihood = math.exp(PL_for + PL_against)

        if n_against > 0:
            self.z_cc_against = (self.z_cc_against / n_against) ** (1 / power)
        if n_for > 0:
            self.z_cc_for = (self.z_cc_for / n_for) ** (1 / power)
        self.z_cc_net = self.z_cc_for - self.z_cc_against
        self.confidence = 0

    def __str__(self):
        """Return a string representation of the subgroup scores.

        Returns:
          str:

        """
        return "%s %.3f %.2f %.2f %.2f %.2f %.2f" % (
            self.subgroup["best_subsym"].space_group_info(),
            self.likelihood,
            self.z_cc_net,
            self.z_cc_for,
            self.z_cc_against,
            self.cc_for.coefficient(),
            self.cc_against.coefficient(),
        )

    def as_dict(self):
        """Return a dictionary representation of the subgroup scoring.

        The dictionary will contain the following keys:
          - patterson_group: The current subgroup
          - likelihood: The likelihood of the subgroup being correct
          - confidence: The confidence of the subgroup being correct
          - z_cc_for: The combined Z-scores for all symmetry elements present in the
            subgroup
          - z_cc_against: The combined Z-scores for all symmetry elements present in
            the lattice group but not in the subgroup
          - z_cc_net: The net Z-score, i.e. z_cc_for - z_cc_against
          - cc_for: The overall correlation coefficient for all symmetry elements
            present in the subgroup
          - cc_against: The overall correlation coefficient for all symmetry
            elements present in the lattice group but not in the subgroup
          - max_angular_difference: The maximum angular difference between the
            symmetrised unit cell and the P1 unit cell.
          - cb_op: The change of basis operation from the input unit cell to the
            'best' unit cell.

        Returns:
          dict:

        """
        return {
            "patterson_group": self.subgroup["best_subsym"]
            .space_group()
            .type()
            .hall_symbol(),
            "likelihood": self.likelihood,
            "confidence": self.confidence,
            "z_cc_net": "% .2f" % self.z_cc_net,
            "z_cc_for": "% .2f" % self.z_cc_for,
            "z_cc_against": "% .2f" % self.z_cc_against,
            "cc_for": "% .2f" % self.cc_for.coefficient(),
            "cc_against": "% .2f" % self.cc_against.coefficient(),
            "max_angular_difference": "%.1f" % self.subgroup["max_angular_difference"],
            "cb_op": "%s" % (self.subgroup["cb_op_inp_best"]),
        }


class CorrelationCoefficientAccumulator(object):
    """Class for incremental computation of correlation coefficients.

    Uses the single-pass formula for Pearson correlation coefficient:
      https://en.wikipedia.org/wiki/Pearson_correlation_coefficient#For_a_sample

    """

    def __init__(self, x=None, y=None):
        """Initialise a CorrelationCoefficientAccumulator object.

        Args:
          x (list): Optional list of `x` values to initialise the accumulator.
          y (list): Optional list of `y` values to initialise the accumulator.

        """
        self._n = 0
        self._sum_x = 0
        self._sum_y = 0
        self._sum_xy = 0
        self._sum_x_sq = 0
        self._sum_y_sq = 0
        if x is not None and y is not None:
            self.accumulate(x, y)

    def accumulate(self, x, y):
        """Accumulate the `x` and `y` values provided.

        Args:
          x (list): The list of `x` values to accumulate.
          y (list): The list of `y` values to accumulate.

        """
        assert x.size() == y.size()
        self._n += x.size()
        self._sum_x += flex.sum(x)
        self._sum_y += flex.sum(y)
        self._sum_xy += flex.sum(x * y)
        self._sum_x_sq += flex.sum(flex.pow2(x))
        self._sum_y_sq += flex.sum(flex.pow2(y))

    def coefficient(self):
        """Calculate the correlation coefficient.

        Returns:
          float: The correlation coefficient.

        """
        if self._n == 0:
            return 0
        return self.numerator() / self.denominator()

    def n(self):
        """Return the number of values contributing to the correlation coefficient.

        Returns:
          n (int)

        """
        return self._n

    def numerator(self):
        r""""Calculate the numerator of the correlation coefficient formula.

        .. math:: n \sum{x y} - \sum{x} \sum{y}

        Returns:
          float: The value of the numerator.

        """
        return self._n * self._sum_xy - self._sum_x * self._sum_y

    def denominator(self):
        r""""Calculate the denominator of the correlation coefficient formula.

        .. math:: \sqrt{n \sum{x^2} - \sum{x}^2} \sqrt{n \sum{y^2} - \sum{y}^2}

        Returns:
          float: The value of the denominator.

        """
        return math.sqrt(self._n * self._sum_x_sq - self._sum_x ** 2) * math.sqrt(
            self._n * self._sum_y_sq - self._sum_y ** 2
        )

    def __iadd__(self, other):
        """Add together two instances of :class:`CorrelationCoefficientAccumulator`.

        Args:
          other (CorrelationCoefficientAccumulator):
            The :class:`CorrelationCoefficientAccumualator` to add to the current object.

        Returns:
          self (CorrelationCoefficientAccumulator): The current object.

        """
        self._n += other._n
        self._sum_x += other._sum_x
        self._sum_y += other._sum_y
        self._sum_xy += other._sum_xy
        self._sum_x_sq += other._sum_x_sq
        self._sum_y_sq += other._sum_y_sq
        return self


def trunccauchy_pdf(x, a, b, loc=0, scale=1):
    """Calculate a truncated Cauchy probability density function.

    Args:
      x (float): The point at which to calculate the PDF.
      a (float): The lower bound of the truncated distribution.
      b (float): The upper bound of the truncated distribution.
      loc (float): The location parameter for the Cauchy distribution.
      scale (float): The scale parameter for the Cauchy distribution.

    Returns:
      float: The value of the probability density function.

    """
    assert b > a
    rv = scipy.stats.cauchy(loc=loc, scale=scale)
    return rv.pdf(x) / (rv.cdf(b) - rv.cdf(a))

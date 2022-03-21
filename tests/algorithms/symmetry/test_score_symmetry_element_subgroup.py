from __future__ import annotations

import pytest

from cctbx import sgtbx
from cctbx.sgtbx.lattice_symmetry import metric_subgroups

from dials.algorithms.symmetry.cosym._generate_test_data import generate_intensities
from dials.algorithms.symmetry.laue_group import (
    ScoreCorrelationCoefficient,
    ScoreSubGroup,
    ScoreSymmetryElement,
)


def test_score_correlation_coefficient():
    cc = 1
    expected_cc = 1
    sigma_cc = 0.1
    score_cc = ScoreCorrelationCoefficient(cc, sigma_cc, expected_cc)
    assert score_cc.p_s_given_cc == pytest.approx(0.9492499653267421)

    cc = 0.5
    score_cc = ScoreCorrelationCoefficient(cc, sigma_cc, expected_cc)
    assert score_cc.p_s_given_cc == pytest.approx(0.19697339347266518)

    cc = 0.5
    expected_cc = 0.6
    sigma_cc = 0.2
    score_cc = ScoreCorrelationCoefficient(cc, sigma_cc, expected_cc)
    assert score_cc.p_s_given_cc == pytest.approx(0.6178439917879021)

    if 0:
        import numpy as np
        from matplotlib import pyplot as plt

        x = np.linspace(-1, 1)
        values = [(0.1, 1), (0.1, 0.8), (0.2, 1), (0.2, 0.8)]
        fig, axes = plt.subplots(nrows=2, ncols=2)
        for ax, (sigma_cc, expected_cc) in zip(axes.flatten(), values):
            ax.set_title(f"E(CC) = {expected_cc:.1f}, sigma(cc) = {sigma_cc:.1f}")
            p_cc_given_s = []
            p_cc_given_not_s = []
            p_s_given_cc = []
            for _ in x:
                score_cc = ScoreCorrelationCoefficient(_, sigma_cc, expected_cc)
                p_cc_given_s.append(score_cc.p_cc_given_s)
                p_cc_given_not_s.append(score_cc.p_cc_given_not_s)
                p_s_given_cc.append(score_cc.p_s_given_cc)
            ax.plot(x, p_cc_given_s, label="p(CC;S)")
            ax.plot(x, p_cc_given_not_s, label="p(CC;!S)")
            ax.plot(x, p_s_given_cc, label="p(S;CC)")
        ylim = max(ax.get_ylim()[1] for ax in axes.flatten())
        for ax in axes.flatten():
            ax.legend()
            ax.set_ylim(ax.get_ylim()[0], ylim)
        plt.show()


@pytest.mark.parametrize("space_group", ["P2", "P3", "P6", "R3:h", "I23"][:])
def test_score_symmetry_element_subgroup(space_group):
    sgi = sgtbx.space_group_info(symbol=space_group)
    cs = sgi.any_compatible_crystal_symmetry(volume=10000)
    cs = cs.best_cell()
    cs = cs.minimum_cell()
    intensities = (
        generate_intensities(cs, d_min=1.0)
        .generate_bijvoet_mates()
        .set_observation_type_xray_intensity()
    )
    intensities = intensities.expand_to_p1()

    subgroups = metric_subgroups(
        intensities.crystal_symmetry(), max_delta=2.0, bravais_types_only=False
    )
    intensities = intensities.change_basis(subgroups.cb_op_inp_minimum)
    cs = cs.change_basis(subgroups.cb_op_inp_minimum)

    cb_op_inp_best = subgroups.result_groups[0]["cb_op_inp_best"]
    lattice_group = subgroups.result_groups[0]["best_subsym"].space_group()
    lattice_group = lattice_group.change_basis(cb_op_inp_best.inverse())

    lattice_group = subgroups.result_groups[0]["subsym"].space_group()

    sym_op_scores = []
    for sym_op in lattice_group.smx():
        if sym_op.r().info().sense() < 0:
            continue
        score = ScoreSymmetryElement(intensities, sym_op, 1.0, 1.0)
        sym_op_scores.append(score)
        if sym_op in cs.space_group():
            assert score.likelihood > 0.9
            assert score.cc.coefficient() > 0.9
        else:
            assert score.likelihood < 0.2
            assert score.cc.coefficient() < 0.3

    subgroup_scores = [
        ScoreSubGroup(subgrp, sym_op_scores) for subgrp in subgroups.result_groups
    ]
    total_likelihood = sum(score.likelihood for score in subgroup_scores)
    for score in subgroup_scores:
        score.likelihood /= total_likelihood
    true_patterson_group = (
        cs.space_group_info()
        .as_reference_setting()
        .group()
        .build_derived_patterson_group()
    )
    for score in subgroup_scores:
        if score.subgroup["best_subsym"].space_group() == true_patterson_group:
            assert score.likelihood > 0.8
        else:
            assert score.likelihood < 0.1

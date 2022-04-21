from __future__ import annotations

import pytest

from cctbx import crystal, miller, sgtbx

from dials.algorithms.symmetry.cosym._generate_test_data import generate_intensities
from dials.algorithms.symmetry.laue_group import LaueGroupAnalysis


def generate_fake_intensities(crystal_symmetry):
    intensities = (
        generate_intensities(crystal_symmetry, d_min=1.0)
        .generate_bijvoet_mates()
        .set_observation_type_xray_intensity()
    )
    intensities = intensities.expand_to_p1()
    # needed to give vaguely sensible E_cc_true values
    intensities = intensities.customized_copy(sigmas=intensities.data() / 50)
    intensities.set_info(miller.array_info(source="fake", source_type="mtz"))
    return intensities


@pytest.mark.parametrize("space_group", ["P2", "P3", "P6", "R3:h", "I23"][:])
def test_determine_space_group(space_group):
    sgi = sgtbx.space_group_info(symbol=space_group)
    sg = sgi.group()
    cs = sgi.any_compatible_crystal_symmetry(volume=10000)
    cs = cs.best_cell()
    cs = cs.minimum_cell()

    intensities = generate_fake_intensities(cs)
    result = LaueGroupAnalysis([intensities], normalisation=None)
    print(result)
    assert (
        result.best_solution.subgroup["best_subsym"].space_group()
        == sg.build_derived_patterson_group()
    )
    assert result.best_solution.likelihood > 0.8
    for score in result.subgroup_scores[1:]:
        assert score.likelihood < 0.1


def test_determine_space_group_best_monoclinic_beta():
    cs = crystal.symmetry(
        unit_cell=(44.66208171, 53.12629403, 62.53397661, 64.86329707, 78.27343894, 90),
        space_group_symbol="C 1 2/m 1 (z,x+y,-2*x)",
    )
    cs = cs.best_cell().primitive_setting()
    intensities = generate_fake_intensities(cs)
    result = LaueGroupAnalysis([intensities], normalisation=None)
    print(result)
    assert result.best_solution.subgroup["best_subsym"].is_similar_symmetry(
        crystal.symmetry(
            unit_cell=(44.66208171, 53.12629403, 111.9989451, 90, 99.89337396, 90),
            space_group_symbol="I 1 2/m 1",
        )
    )
    d = result.as_dict()
    assert cs.change_basis(
        sgtbx.change_of_basis_op(d["subgroup_scores"][0]["cb_op"])
    ).is_similar_symmetry(result.best_solution.subgroup["best_subsym"])
    result = LaueGroupAnalysis(
        [intensities], best_monoclinic_beta=False, normalisation=None
    )
    print(result)
    assert result.best_solution.subgroup["best_subsym"].is_similar_symmetry(
        crystal.symmetry(
            unit_cell=(113.2236274, 53.12629403, 44.66208171, 90, 102.9736126, 90),
            space_group_symbol="C 1 2/m 1",
        )
    )
    d = result.as_dict()
    assert cs.change_basis(
        sgtbx.change_of_basis_op(d["subgroup_scores"][0]["cb_op"])
    ).is_similar_symmetry(result.best_solution.subgroup["best_subsym"])

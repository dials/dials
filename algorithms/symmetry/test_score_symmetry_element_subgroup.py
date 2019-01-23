from __future__ import absolute_import, division, print_function

import pytest

from cctbx import sgtbx
from cctbx.sgtbx.lattice_symmetry import metric_subgroups
from cctbx.sgtbx.subgroups import subgroups
from dials.algorithms.symmetry.cosym._generate_test_data import generate_intensities
from dials.algorithms.symmetry.determine_space_group import ScoreSymmetryElement
from dials.algorithms.symmetry.determine_space_group import ScoreSubGroup

@pytest.mark.parametrize('space_group', ['P2', 'P3', 'P6', 'R3:h', 'I23'][:])
def test_score_symmetry_element_subgroup(space_group):
  sgi = sgtbx.space_group_info(symbol=space_group)
  sg = sgi.group()
  cs = sgi.any_compatible_crystal_symmetry(volume=10000)
  cs = cs.best_cell()
  cs = cs.minimum_cell()
  intensities = generate_intensities(
    cs, d_min=1.0).generate_bijvoet_mates().set_observation_type_xray_intensity()
  intensities = intensities.expand_to_p1()

  subgroups = metric_subgroups(
    intensities.crystal_symmetry(), max_delta=2.0, bravais_types_only=False)
  cb_op_inp_best = subgroups.result_groups[0]['cb_op_inp_best']
  lattice_group = subgroups.result_groups[0]['best_subsym'].space_group()
  lattice_group = lattice_group.change_basis(cb_op_inp_best.inverse())

  sym_op_scores = []
  for sym_op in lattice_group.smx():
    if sym_op.r().info().sense() < 0: continue
    score = ScoreSymmetryElement(intensities, sym_op, 1.0, 1.0)
    sym_op_scores.append(score)
    if sym_op in cs.space_group():
      assert score.likelihood > 0.9
      assert score.cc.coefficient() > 0.9
    else:
      assert score.likelihood < 0.2
      assert score.cc.coefficient() < 0.3

  subgroup_scores = [
    ScoreSubGroup(subgrp, sym_op_scores)
    for subgrp in subgroups.result_groups]
  total_likelihood = sum(score.likelihood for score in subgroup_scores)
  for score in subgroup_scores: score.likelihood /= total_likelihood
  true_patterson_group = cs.space_group_info().as_reference_setting().group().build_derived_patterson_group()
  for score in subgroup_scores:
    if score.subgroup['best_subsym'].space_group() == true_patterson_group:
      assert score.likelihood > 0.8
    else:
      assert score.likelihood < 0.1

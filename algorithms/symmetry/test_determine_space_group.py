from __future__ import absolute_import, division, print_function

import pytest

from cctbx import miller
from cctbx import sgtbx
from dials.algorithms.symmetry.cosym.generate_test_data import generate_intensities
from dials.algorithms.symmetry.determine_space_group import determine_space_group

@pytest.mark.parametrize('space_group', ['P2', 'P3', 'P6', 'R3:h', 'I23'][:])
def test_determine_space_group(space_group):

  sgi = sgtbx.space_group_info(symbol=space_group)
  sg = sgi.group()
  cs = sgi.any_compatible_crystal_symmetry(volume=10000)
  cs = cs.best_cell()
  cs = cs.minimum_cell()
  intensities = generate_intensities(
    cs, d_min=1.0).generate_bijvoet_mates().set_observation_type_xray_intensity()
  intensities = intensities.expand_to_p1()
  # needed to give vaguely sensible E_cc_true values
  intensities = intensities.customized_copy(sigmas=intensities.data()/50)
  intensities.set_info(miller.array_info(
    source='fake',
    source_type='mtz'
  ))
  result = determine_space_group([intensities], normalisation=None)
  #import logging
  #logging.basicConfig(level=logging.INFO)
  #result.show()
  assert result.best_solution.subgroup['best_subsym'].space_group() == sg.build_derived_patterson_group()
  assert result.best_solution.likelihood > 0.8
  for score in result.subgroup_scores[1:]:
    assert score.likelihood < 0.1

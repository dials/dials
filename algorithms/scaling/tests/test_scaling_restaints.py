from dials.array_family import flex
from dials.algorithms.scaling.scaling_restraints import ScalingRestraints,\
  MultiScalingRestraints
from dials.algorithms.scaling.active_parameter_managers import \
  active_parameter_manager, multi_active_parameter_manager
from dials.algorithms.scaling.model.components.scale_components import \
  SHScaleComponent, SingleScaleFactor

def test_ScalingRestraints():
  """Test for the single scaling restraints manager."""
  abs_comp = SHScaleComponent(flex.double([0.8, 1.0, 0.8]))
  abs_comp.parameter_restraints = flex.double([0.1, 0.1, 0.1])
  scale_comp = SingleScaleFactor(flex.double([1.0]))
  components = {'scale' : scale_comp, 'absorption' : abs_comp}

  apm = active_parameter_manager(components, ['scale', 'absorption'])

  restraints = ScalingRestraints(apm).calculate_restraints()
  abs_restraints = components['absorption'].calculate_restraints()
  assert list(restraints[0]) == list(abs_restraints[0])
  assert list(restraints[1]) == list(abs_restraints[1]) + [0.0]

  jacobian_restraints = ScalingRestraints(apm).calculate_jacobian_restraints()
  abs_restraints = components['absorption'].calculate_jacobian_restraints()

  assert list(jacobian_restraints[0]) == list(abs_restraints[0])
  assert jacobian_restraints[1].n_cols == 4
  assert jacobian_restraints[1].n_rows == 3
  for i in range(3):
    for j in range(3):
      assert jacobian_restraints[1][i, j] == abs_restraints[1][i, j]
    assert jacobian_restraints[1][i, 3] == 0.0

def test_MultiScalingRestraints():
  """Test for the single scaling restraints manager."""

  # Set up some scaling model compoents and a multi active parameter manager.
  abs_comp = SHScaleComponent(flex.double([0.8, 1.0, 0.8]))
  abs_comp.parameter_restraints = flex.double([0.1, 0.1, 0.1])
  scale_comp = SingleScaleFactor(flex.double([1.0]))
  components_1 = {'scale' : scale_comp, 'absorption' : abs_comp}
  abs_comp_2 = SHScaleComponent(flex.double([0.9, 1.0, 0.9]))
  abs_comp_2.parameter_restraints = flex.double([0.2, 0.2, 0.2])
  scale_comp_2 = SingleScaleFactor(flex.double([1.0]))
  components_2 = {'scale' : scale_comp_2, 'absorption' : abs_comp_2}

  apm = multi_active_parameter_manager([components_1, components_2], [
    ['scale', 'absorption'], ['scale', 'absorption']],
    apm_class=active_parameter_manager)

  # Initialise and test that the restaints are correctly composed from the
  # individual component restraints.
  restraints = MultiScalingRestraints(apm).calculate_restraints()
  abs_restraints_1 = components_1['absorption'].calculate_restraints()
  abs_restraints_2 = components_2['absorption'].calculate_restraints()
  assert list(restraints[0]) == (list(abs_restraints_1[0]) +
    list(abs_restraints_2[0]))
  assert list(restraints[1]) == (list(abs_restraints_1[1]) + [0.0] +
    list(abs_restraints_2[1]) + [0.0])

  # Expected to create an n_abs_params by n_total_params matrix, where the abs
  # jacobian restraints matrices are correctly located with respect to their
  # order in the apm and each on their own row.
  jacobian_restraints = MultiScalingRestraints(apm).calculate_jacobian_restraints()
  abs_restraints_1 = components_1['absorption'].calculate_jacobian_restraints()
  abs_restraints_2 = components_2['absorption'].calculate_jacobian_restraints()

  assert list(jacobian_restraints[0]) == (list(abs_restraints_1[0]) +
    list(abs_restraints_2[0]))
  assert jacobian_restraints[1].n_cols == 8
  assert jacobian_restraints[1].n_rows == 6
  for i in range(3):
    for j in range(3):
      assert jacobian_restraints[1][i, j] == abs_restraints_1[1][i, j]
    assert jacobian_restraints[1][i, 3] == 0.0
  for i in range(3):
    for j in range(3):
      assert jacobian_restraints[1][i+3, j+4] == abs_restraints_2[1][i, j]
    assert jacobian_restraints[1][i+3, 7] == 0.0

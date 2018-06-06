"""
Test for the basis function module.
"""
import pytest
from dials.array_family import flex
from dials.algorithms.scaling.model.components.scale_components import \
  SingleBScaleFactor, SingleScaleFactor
from dials.algorithms.scaling.basis_functions import basis_function
from dials.algorithms.scaling.parameter_handler import \
  scaling_active_parameter_manager

@pytest.fixture
def small_reflection_table():
  """Generate reflection table to test the basis function."""
  reflections = flex.reflection_table()
  reflections['d'] = flex.double([2.0, 0.8, 2.0]) #don't change
  reflections['xyzobs.px.value'] = flex.vec3_double([(0.0, 0.0, 0.0),
    (0.0, 0.0, 5.0), (0.0, 0.0, 10.0)])
  return reflections

def test_basis_function(small_reflection_table):
  """Test for the basis function class. This calculates scale factors and
  derivatives for reflections based on the model components."""

  # To test the basis function, need a scaling active parameter manager - to set
  # this up we need a components dictionary with some reflection data.

  # Let's use KB model components for simplicity.
  rt = small_reflection_table
  components = {'scale' : SingleScaleFactor(flex.double([1.0])), 'decay':
    SingleBScaleFactor(flex.double([0.0]))} #Create empty components.
  for component in components.itervalues():
    component.update_reflection_data(rt) #Add some data to components.

  apm = scaling_active_parameter_manager(components, ['decay', 'scale'])

  # First test that scale factors can be successfully updated.
  # Manually change the parameters in the apm.
  decay = components['decay'] # Define alias
  scale = components['scale'] # Define alias
  # Note, order of params in apm.x depends on order in scaling model components.
  new_B = 1.0
  new_S = 2.0
  apm.set_param_vals(flex.double([new_S, new_B]))
  basis_fn = basis_function(curvatures=True)
  basis_fn.update_scale_factors(apm)

  # Now test that the inverse scale factor is correctly calculated.
  calculated_sfs = basis_fn.calculate_scale_factors(apm)[0]
  assert list(calculated_sfs) == list(new_S * flex.exp(new_B/
    (2.0*(decay.d_values[0]**2))))

  # Now check that the derivative matrix is correctly calculated.
  calc_derivs = basis_fn.calculate_derivatives(apm)[0]
  assert calc_derivs[0, 0] == scale.derivatives[0][0, 0] * decay.inverse_scales[0][0]
  assert calc_derivs[1, 0] == scale.derivatives[0][1, 0] * decay.inverse_scales[0][1]
  assert calc_derivs[2, 0] == scale.derivatives[0][2, 0] * decay.inverse_scales[0][2]
  assert calc_derivs[0, 1] == decay.derivatives[0][0, 0] * scale.inverse_scales[0][0]
  assert calc_derivs[1, 1] == decay.derivatives[0][1, 0] * scale.inverse_scales[0][1]
  assert calc_derivs[2, 1] == decay.derivatives[0][2, 0] * scale.inverse_scales[0][2]

  # Test that the curvatures matrix is correctly composed.
  calc_curvs = basis_fn.calculate_curvatures(apm)[0]
  assert calc_curvs[0, 0] == scale.curvatures[0][0, 0] * decay.inverse_scales[0][0]
  assert calc_curvs[1, 0] == scale.curvatures[0][1, 0] * decay.inverse_scales[0][1]
  assert calc_curvs[2, 0] == scale.curvatures[0][2, 0] * decay.inverse_scales[0][2]
  assert calc_curvs[0, 1] == decay.curvatures[0][0, 0] * scale.inverse_scales[0][0]
  assert calc_curvs[1, 1] == decay.curvatures[0][1, 0] * scale.inverse_scales[0][1]
  assert calc_curvs[2, 1] == decay.curvatures[0][2, 0] * scale.inverse_scales[0][2]

  # Repeat the test when there is only one active parameter.
  # First reset the parameters
  components['decay'].parameters = flex.double([0.0])
  components['scale'].parameters = flex.double([1.0])
  components['decay'].calculate_scales_and_derivatives()
  components['scale'].calculate_scales_and_derivatives()

  # Now generate a parameter manager for a single component.
  apm = scaling_active_parameter_manager(components, ['scale'])
  new_S = 2.0
  apm.set_param_vals(flex.double(components['scale'].n_params, new_S))
  basis_fn = basis_function().calculate_scales_and_derivatives(apm)
  #basis_fn = basis_func.calculate_scales_and_derivatives() # All in one alternative call.

  # Test that the scales and derivatives were correctly calculated
  assert list(basis_fn[0][0]) == list([new_S] *
    components['scale'].inverse_scales[0].size())
  assert basis_fn[1][0][0, 0] == components['scale'].derivatives[0][0, 0]
  assert basis_fn[1][0][1, 0] == components['scale'].derivatives[0][1, 0]
  assert basis_fn[1][0][2, 0] == components['scale'].derivatives[0][2, 0]

  apm = scaling_active_parameter_manager(components, [])
  basis_fn = basis_function(curvatures=True)
  _, d, c = basis_fn.calculate_scales_and_derivatives(apm)
  assert d is None
  assert c is None

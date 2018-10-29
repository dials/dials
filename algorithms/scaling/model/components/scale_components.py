""""
Classes that each define one component of a scaling model.

These classes hold the component parameters and a component of
the inverse scale factors, as well as methods to calculate and
set the inverse scale factors and derivatives.
"""
import abc
from dials.array_family import flex
from scitbx import sparse


class ScaleComponentBase(object):
  """Base scale component class.

  This defines an interface to access the parameters, the component
  of the inverse scale factor and it's derivatives with respect to
  the parameters. Scale components derived from the base class are
  designed to be instantiated by a ScalingModel class, by supplying
  an initial array of parameters and optionally the current estimated
  standard deviations. The relevant data from a reflection table is
  added later by a Scaler using the update_reflection_data method.
  This behaviour allows data to easily be added/changed after selecting
  subsets of the data."""

  __metaclass__ = abc.ABCMeta

  def __init__(self, initial_values, parameter_esds=None):
    self._parameters = initial_values
    self._parameter_esds = parameter_esds
    self._n_params = len(self._parameters)
    self._var_cov = None
    # Store the below properties as lists, to allow splitting into blocks
    # for multiprocessing.
    self._n_refl = []
    self._inverse_scales = []
    self._derivatives = []
    self._curvatures = []
    self._parameter_restraints = None

  @property
  def parameter_restraints(self):
    """Restraint weights for the component parameters."""
    return self._parameter_restraints

  @parameter_restraints.setter
  def parameter_restraints(self, restraints):
    assert restraints.size() == self.parameters.size()
    """Set Restraint weights for the component parameters."""
    self._parameter_restraints = restraints

  @property
  def n_params(self):
    """The number of parameters that parameterise the component."""
    return self._n_params

  @property
  def parameters(self):
    """The parameters of the component."""
    return self._parameters

  @parameters.setter
  def parameters(self, values):
    if len(values) != len(self._parameters):
      assert 0, '''attempting to set a new set of parameters of different
      length than previous assignment: was %s, attempting %s''' % (
        len(self._parameters), len(values))
    self._parameters = values

  @property
  def parameter_esds(self):
    """The estimated standard deviations of the parameters."""
    return self._parameter_esds

  @parameter_esds.setter
  def parameter_esds(self, esds):
    assert len(esds) == len(self._parameters)
    self._parameter_esds = esds

  def calculate_restraints(self):
    """Calculate residual and gradient restraints for the component."""
    return None

  def calculate_jacobian_restraints(self):
    """Calculate residual and jacobian restraints for the component."""
    return None

  @property
  def var_cov_matrix(self):
    """The variance/covariance matrix of the parameters."""
    return self._var_cov

  @var_cov_matrix.setter
  def var_cov_matrix(self, var_cov):
    self._var_cov = var_cov

  @property
  def n_refl(self):
    """The current number of reflections associated with this component."""
    return self._n_refl

  @property
  def inverse_scales(self):
    """The inverse scale factors associated with the data."""
    return self._inverse_scales

  @inverse_scales.setter
  def inverse_scales(self, new_inverse_scales):
    self._inverse_scales = new_inverse_scales

  @property
  def derivatives(self):
    """A spares matrix of the derivatives of the inverse scale
    factors with respect to the component parameters."""
    return self._derivatives

  @property
  def curvatures(self):
    """A spares matrix of the curvatures of the inverse scale
    factors with respect to the component parameters."""
    return self._curvatures

  @abc.abstractmethod
  def update_reflection_data(self, reflection_table, selection=None,
    block_selections=None):
    """Add or change the relevant reflection data for the component.

    No restrictions should be placed on the data remaining the same
    size, to allow easy changing of the data contained during scaling.
    The input data will be specific to the component."""

  @abc.abstractmethod
  def calculate_scales_and_derivatives(self, block_id=0, curvatures=False):
    """Use the component parameters to calculate and set
    self._inverse_scales and self._derivatives."""


class SingleScaleFactor(ScaleComponentBase):
  """A model component consisting of a single global scale parameter.

  The inverse scale factor for every reflection is the parameter
  value and the derivatives are all 1.0."""

  def __init__(self, initial_values, parameter_esds=None):
    assert len(initial_values) == 1, '''This model component is only designed
      for a single global scale component'''
    super(SingleScaleFactor, self).__init__(initial_values, parameter_esds)

  def update_reflection_data(self, reflection_table, selection=None,
    block_selections=None):
    """Add reflection data to the component, only n_reflections needed."""
    self._n_refl = []
    if selection:
      reflections = reflection_table.select(selection)
    else:
      reflections = reflection_table
    if block_selections:
      block_selection_list = block_selections
      for sel in block_selection_list:
        self._n_refl.append(reflections.select(sel).size())
    else:
      self._n_refl.append(reflections.size())

  def calculate_scales_and_derivatives(self, block_id=0, curvatures=False):
    scales = flex.double(self.n_refl[block_id],
        self._parameters[0])
    derivatives = sparse.matrix(self.n_refl[block_id], 1)
    for i in range(self.n_refl[block_id]):
      derivatives[i, 0] = 1.0
    if curvatures:
      curvatures = sparse.matrix(self.n_refl[block_id], 1)#curvatures are all zero.
      return scales, derivatives, curvatures
    return scales, derivatives

class SingleBScaleFactor(ScaleComponentBase):
  """A model component for a single global B-factor parameter.

  The inverse scale factor for each reflection is given by
  S = exp(B/(2 * d^2)), the derivatives are S/(2 * d^2)."""

  def __init__(self, initial_values, parameter_esds=None):
    super(SingleBScaleFactor, self).__init__(initial_values, parameter_esds)
    self._d_values = []

  @property
  def d_values(self):
    """The current set of d-values associated with this component."""
    return self._d_values

  def update_reflection_data(self, reflection_table, selection=None,
    block_selections=None):
    """"Add reflection data to the component, only the d-values and number
    of reflections are needed."""
    self._n_refl = []
    self._d_values = []
    if selection:
      reflections = reflection_table.select(selection)
    else:
      reflections = reflection_table
    if block_selections:
      block_selection_list = block_selections
      for i, sel in enumerate(block_selection_list):
        self._d_values.append(reflections['d'].select(sel))
        self._n_refl.append(self._d_values[i].size())
    else:
      self._d_values.append(reflections['d'])
      self._n_refl.append(reflections.size())

  def calculate_scales_and_derivatives(self, block_id=0, curvatures=False):
    scales = flex.exp(flex.double(
      [self._parameters[0]] * self._n_refl[block_id])
      / (2.0 * (self._d_values[block_id]*self._d_values[block_id])))
    derivatives = sparse.matrix(self._n_refl[block_id], 1)
    for i in range(self._n_refl[block_id]):
      derivatives[i, 0] = (scales[i]
        / (2.0 * (self._d_values[block_id][i]*self._d_values[block_id][i])))
    if curvatures:
      curvatures = sparse.matrix(self.n_refl[block_id], 1) #curvatures are all zero.
      for i in range(self._n_refl[block_id]):
        curvatures[i, 0] = (scales[i]
          / ((2.0 * (self._d_values[block_id][i]**2))**2))
      return scales, derivatives, curvatures
    return scales, derivatives

class SHScaleComponent(ScaleComponentBase):
  """A model component for a spherical harmonic absorption correction.

  This component uses a set of spherical harmonic functions to define
  an absorption surface for the crystal. A matrix of spherical harmonic
  coefficients for the data is stored in self._harmonic_values and is
  used to calculate the scales and derivatives.
  The scale is given by S = 1 + (sum_l sum_m Clm * Ylm) where Clm are
  the model parameters and Ylm are the spherical harmonic coefficients,
  the derivatives are then simply the coefficients Ylm."""

  def __init__(self, initial_values, parameter_esds=None):
    super(SHScaleComponent, self).__init__(initial_values, parameter_esds)
    self._harmonic_values = []
    self._sph_harm_table = None
    self._parameter_restraints = None

  @property
  def harmonic_values(self):
    """A matrix of harmonic coefficients for the data."""
    return self._harmonic_values

  @property
  def sph_harm_table(self):
    """A matrix of the full harmonic coefficient for a reflection table."""
    return self._sph_harm_table

  @sph_harm_table.setter
  def sph_harm_table(self, sht):
    """Set the spherical harmonic table."""
    self._sph_harm_table = sht

  def calculate_restraints(self):
    residual = self.parameter_restraints * (self._parameters*self._parameters)
    gradient = 2.0 * self.parameter_restraints * self._parameters
    return residual, gradient

  def calculate_jacobian_restraints(self):
    jacobian = sparse.matrix(self.n_params, self.n_params)
    for i in range(self.n_params):
      jacobian[i, i] = +1.0
    return self._parameters, jacobian, self._parameter_restraints

  def update_reflection_data(self, _, selection=None, block_selections=None):
    """Update the spherical harmonic coefficients."""
    self._n_refl = []
    self._harmonic_values = []
    if selection:
      sel_sph_harm_table = self.sph_harm_table.select_columns(
        selection.iselection())
      if block_selections:
        block_selection_list = block_selections
        for i, sel in enumerate(block_selection_list):
          block_sph_harm_table = sel_sph_harm_table.select_columns(sel)
          self._harmonic_values.append(block_sph_harm_table.transpose())
          self.n_refl.append(self._harmonic_values[i].n_rows)
      else:
        self._harmonic_values.append(sel_sph_harm_table.transpose())
        self._n_refl.append(self._harmonic_values[0].n_rows)
    else:
      self._harmonic_values.append(self.sph_harm_table.transpose())
      self._n_refl.append(self._harmonic_values[0].n_rows)

  def calculate_scales_and_derivatives(self, block_id=0, curvatures=False):
    abs_scale = flex.double(self._harmonic_values[block_id].n_rows, 1.0
      ) #Unity term
    for i, col in enumerate(self._harmonic_values[block_id].cols()):
      abs_scale += flex.double(col.as_dense_vector() * self._parameters[i])
    if curvatures:
      curvatures = sparse.matrix(
        abs_scale.size(), self._n_params)
      return abs_scale, self._harmonic_values[block_id], curvatures
    return abs_scale, self._harmonic_values[block_id]

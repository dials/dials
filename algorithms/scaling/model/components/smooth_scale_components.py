""""
Classes that each define a smoothly varying component of a scaling model.

These classes use a gaussian smoother (1D, 2D or 3D) to calculate the
inverse scale factors and derivatives with respect to the component
parameters.
"""
from __future__ import absolute_import, division, print_function
from math import floor, ceil
from scitbx import sparse
from dials.array_family import flex
from dials.algorithms.scaling.model.components.scale_components import (
    ScaleComponentBase,
)
from dials_scaling_ext import row_multiply
from dials_refinement_helpers_ext import GaussianSmoother as GS1D
from dials_refinement_helpers_ext import GaussianSmoother2D as GS2D
from dials_refinement_helpers_ext import GaussianSmoother3D as GS3D


# The following gaussian smoother classes make the implementation
# consistent with that used in dials.refinement.


class GaussianSmoother1D(GS1D):
    """A 1D Gaussian smoother."""

    def value_weight(self, x, value):
        """Return the value, weight and sumweight at a single point."""
        result = super(GaussianSmoother1D, self).value_weight(x, value)
        return (result.get_value(), result.get_weight(), result.get_sumweight())

    def multi_value_weight(self, x, value):
        """Return the value, weight and sumweight at multiple points."""
        result = super(GaussianSmoother1D, self).multi_value_weight(x, value)
        return (result.get_value(), result.get_weight(), result.get_sumweight())

    def positions(self):
        """Return the smoother positions."""
        return list(super(GaussianSmoother1D, self).positions())


class GaussianSmoother2D(GS2D):
    """A 2D Gaussian smoother."""

    def value_weight(self, x, y, value):
        """Return the value, weight and sumweight at a single point."""
        result = super(GaussianSmoother2D, self).value_weight(x, y, value)
        return (result.get_value(), result.get_weight(), result.get_sumweight())

    def multi_value_weight(self, x, y, value):
        """Return the value, weight and sumweight at multiple points."""
        result = super(GaussianSmoother2D, self).multi_value_weight(x, y, value)
        return (result.get_value(), result.get_weight(), result.get_sumweight())

    def x_positions(self):
        """Return the smoother x-positions."""
        return list(super(GaussianSmoother2D, self).x_positions())

    def y_positions(self):
        """Return the smoother y-positions."""
        return list(super(GaussianSmoother2D, self).y_positions())


class GaussianSmoother3D(GS3D):
    """A 3D Gaussian smoother."""

    def value_weight(self, x, y, z, value):
        """Return the value, weight and sumweight at a single point."""
        result = super(GaussianSmoother3D, self).value_weight(x, y, z, value)
        return (result.get_value(), result.get_weight(), result.get_sumweight())

    def multi_value_weight(self, x, y, z, value):
        """Return the value, weight and sumweight at multiple points."""
        result = super(GaussianSmoother3D, self).multi_value_weight(x, y, z, value)
        return (result.get_value(), result.get_weight(), result.get_sumweight())

    def x_positions(self):
        """Return the smoother x-positions."""
        return list(super(GaussianSmoother3D, self).x_positions())

    def y_positions(self):
        """Return the smoother y-positions."""
        return list(super(GaussianSmoother3D, self).y_positions())

    def z_positions(self):
        """Return the smoother z-positions."""
        return list(super(GaussianSmoother3D, self).z_positions())


class SmoothMixin(object):
    """Mixin class for smooth scale factor components.

    This uses a Gaussian smoother to calculate scales and derivatives
    based on the parameters and a have a set of normalised_values
    associated with the data."""

    def __init__(self):
        self._Vr = 1.0
        self._smoother = None

    @property
    def value(self):
        """Extra access to the parameters for the gaussian smoother."""
        return self._parameters

    @property
    def smoother(self):
        """The Gaussian smoother."""
        return self._smoother

    @staticmethod
    def nparam_to_val(n_params):
        """Convert the number of parameters to the required input value
        for the smoother."""
        assert (
            n_params >= 2
        ), """cannot initialise a smooth scale factor
      for <2 parameters."""
        if n_params == 2 or n_params == 3:
            return n_params - 1
        return n_params - 2


class SmoothScaleComponent1D(ScaleComponentBase, SmoothMixin):
    """A smoothly varying scale component in one dimension."""

    null_parameter_value = 1.0

    def __init__(self, initial_values, parameter_esds=None):
        super(SmoothScaleComponent1D, self).__init__(initial_values, parameter_esds)
        self._normalised_values = []

    def set_new_parameters(self, new_parameters):
        """Set new parameters of a different length i.e. after batch handling"""
        self._parameters = new_parameters
        self._parameter_esds = None
        self._n_params = len(self._parameters)

    @property
    def normalised_values(self):
        """This is a list of the relevant data needed to calculate the
        inverse scale factors, normalised to give 'normalised coordinate
        values' that fit in the range of the smoother parameters, which
        are defined as a 1D array at normalised coordinates separated by
        a spacing of 1."""
        return self._normalised_values

    @ScaleComponentBase.data.setter
    def data(self, data):
        assert set(data.keys()) == {"x"}, set(data.keys())
        self._data = data

    def update_reflection_data(self, selection=None, block_selections=None):
        """Set the normalised coordinate values and configure the smoother."""
        self._normalised_values = []
        self._n_refl = []
        normalised_values = self.data["x"]
        if selection:
            normalised_values = normalised_values.select(selection)
        # Make sure zeroed correctly.
        normalised_values = normalised_values - flex.min(normalised_values)
        phi_range_deg = [
            floor(round(flex.min(normalised_values), 10)),
            ceil(round(flex.max(normalised_values), 10)),
        ]
        self._smoother = GaussianSmoother1D(
            phi_range_deg, self.nparam_to_val(self._n_params)
        )
        if block_selections:
            block_selection_list = block_selections
            for i, sel in enumerate(block_selection_list):
                self._normalised_values.append(normalised_values.select(sel))
                self._n_refl.append(self._normalised_values[i].size())
        else:
            self._normalised_values.append(normalised_values)
            self._n_refl.append(normalised_values.size())

    def calculate_scales_and_derivatives(self, block_id=0):
        if self._n_refl[block_id] > 1:
            value, weight, sumweight = self._smoother.multi_value_weight(
                self._normalised_values[block_id], self.value
            )
            inv_sw = 1.0 / sumweight
            dv_dp = row_multiply(weight, inv_sw)
        elif self._n_refl[block_id] == 1:
            value, weight, sumweight = self._smoother.value_weight(
                self._normalised_values[block_id][0], self.value
            )
            dv_dp = sparse.matrix(1, weight.size)
            b = flex.double(weight.as_dense_vector() / sumweight)
            b.reshape(flex.grid(1, b.size()))
            dv_dp.assign_block(b, 0, 0)
            value = flex.double(1, value)
        else:
            return flex.double([]), sparse.matrix(0, 0)
        return value, dv_dp

    def calculate_scales(self, block_id=0):
        """"Only calculate the scales if needed, for performance."""
        if self._n_refl[block_id] > 1:
            value, _, __ = self._smoother.multi_value_weight(
                self._normalised_values[block_id], self.value
            )
        elif self._n_refl[block_id] == 1:
            value, _, __ = self._smoother.value_weight(
                self._normalised_values[block_id][0], self.value
            )
            value = flex.double(1, value)
        else:
            value = flex.double([])
        return value


class SmoothBScaleComponent1D(SmoothScaleComponent1D):
    """Subclass of SmoothScaleComponent1D to implement a smoothly
    varying B-factor correction."""

    null_parameter_value = 0.0

    def __init__(self, initial_values, parameter_esds=None):
        super(SmoothBScaleComponent1D, self).__init__(initial_values, parameter_esds)
        self._d_values = []

    @property
    def d_values(self):
        """The current set of d-values associated with this component."""
        return self._d_values

    @ScaleComponentBase.data.setter
    def data(self, data):
        assert set(data.keys()) == {"x", "d"}, set(data.keys())
        self._data = data

    def update_reflection_data(self, selection=None, block_selections=None):
        super(SmoothBScaleComponent1D, self).update_reflection_data(
            selection, block_selections
        )
        self._d_values = []
        data = self.data["d"]
        if selection:
            data = data.select(selection)
        if block_selections:
            for sel in block_selections:
                self._d_values.append(data.select(sel))
        else:
            self._d_values.append(data)

    def calculate_scales_and_derivatives(self, block_id=0):
        sdctuple = super(
            SmoothBScaleComponent1D, self
        ).calculate_scales_and_derivatives(block_id)
        if self._n_refl[block_id] == 0:
            return flex.double([]), sparse.matrix(0, 0)
        prefac = 1.0 / (2.0 * (self._d_values[block_id] * self._d_values[block_id]))
        s = flex.exp(sdctuple[0] * prefac)
        d = row_multiply(sdctuple[1], s * prefac)
        return s, d

    def calculate_scales(self, block_id=0):
        s = super(SmoothBScaleComponent1D, self).calculate_scales(block_id)
        return flex.exp(s / (2.0 * (self._d_values[block_id] ** 2)))

    def calculate_restraints(self):
        residual = self.parameter_restraints * (self._parameters * self._parameters)
        gradient = 2.0 * self.parameter_restraints * self._parameters
        return residual, gradient

    def calculate_jacobian_restraints(self):
        jacobian = sparse.matrix(self.n_params, self.n_params)
        for i in range(self.n_params):
            jacobian[i, i] = +1.0
        return self._parameters, jacobian, self.parameter_restraints


class SmoothScaleComponent2D(ScaleComponentBase, SmoothMixin):
    """Implementation of a 2D array-based smoothly varying scale factor.

    A 2d array of parameters is defined, and the scale factor at fractional
    coordinates is calculated as smoothly varying based on the distance to
    the nearby parameters as calculated in the GaussianSmoother2D. The
    initial values are passed as a 1D array, and shape is a 2-tuple
    indicating the number of parameters in each dimension."""

    null_parameter_value = 1.0

    def __init__(self, initial_values, shape, parameter_esds=None):
        assert len(initial_values) == (
            shape[0] * shape[1]
        ), """The shape
    information to initialise a 2D smoother is inconsistent with the length
    of the initial parameter list."""
        super(SmoothScaleComponent2D, self).__init__(initial_values, parameter_esds)
        self._n_x_params = shape[0]
        self._n_y_params = shape[1]
        self._normalised_x_values = None
        self._normalised_y_values = None

    @ScaleComponentBase.data.setter
    def data(self, data):
        assert set(data.keys()) == {"x", "y"}, set(data.keys())
        self._data = data

    def set_new_parameters(self, new_parameters, shape):
        """Set new parameters of a different length i.e. after batch handling"""
        assert len(new_parameters) == shape[0] * shape[1]
        self._parameters = new_parameters
        self._parameter_esds = None
        self._n_params = len(self._parameters)
        self._n_x_params = shape[0]
        self._n_y_params = shape[1]

    @property
    def n_x_params(self):
        """The number of parameters that parameterise the x-component."""
        return self._n_x_params

    @property
    def n_y_params(self):
        """The number of parameters that parameterise the y-component."""
        return self._n_y_params

    @property
    def normalised_x_values(self):
        """The normalised coordinate values in the first dimension."""
        return self._normalised_x_values

    @property
    def normalised_y_values(self):
        """The normalised coordinate values in the second dimension."""
        return self._normalised_y_values

    def update_reflection_data(self, selection=None, block_selections=None):
        """control access to setting all of reflection data at once"""

        self._normalised_x_values = []
        self._normalised_y_values = []
        self._n_refl = []
        normalised_x_values = self.data["x"]
        normalised_y_values = self.data["y"]
        if selection:
            normalised_x_values = normalised_x_values.select(selection)
            normalised_y_values = normalised_y_values.select(selection)
        normalised_x_values = normalised_x_values - flex.min(normalised_x_values)
        normalised_y_values = normalised_y_values - flex.min(normalised_y_values)
        x_range = [
            floor(round(flex.min(normalised_x_values), 10)),
            ceil(round(flex.max(normalised_x_values), 10)),
        ]
        y_range = [
            floor(round(flex.min(normalised_y_values), 10)),
            ceil(round(flex.max(normalised_y_values), 10)),
        ]
        self._smoother = GaussianSmoother2D(
            x_range,
            self.nparam_to_val(self._n_x_params),
            y_range,
            self.nparam_to_val(self._n_y_params),
        )
        if block_selections:
            for i, sel in enumerate(block_selections):
                self._normalised_x_values.append(normalised_x_values.select(sel))
                self._normalised_y_values.append(normalised_y_values.select(sel))
                self._n_refl.append(self._normalised_x_values[i].size())
        else:
            self._normalised_x_values.append(normalised_x_values)
            self._normalised_y_values.append(normalised_y_values)
            self._n_refl.append(normalised_x_values.size())

    def calculate_scales_and_derivatives(self, block_id=0):
        if self._n_refl[block_id] > 1:
            value, weight, sumweight = self._smoother.multi_value_weight(
                self._normalised_x_values[block_id],
                self._normalised_y_values[block_id],
                self.value,
            )
            inv_sw = 1.0 / sumweight
            dv_dp = row_multiply(weight, inv_sw)
        elif self._n_refl[block_id] == 1:
            value, weight, sumweight = self._smoother.value_weight(
                self._normalised_x_values[block_id][0],
                self._normalised_y_values[block_id][0],
                self.value,
            )
            dv_dp = sparse.matrix(1, weight.size)
            b = flex.double(weight.as_dense_vector() / sumweight)
            b.reshape(flex.grid(1, b.size()))
            dv_dp.assign_block(b, 0, 0)
            value = flex.double(1, value)
        else:
            return flex.double([]), sparse.matrix(0, 0)
        return value, dv_dp

    def calculate_scales(self, block_id=0):
        """Only calculate the scales if needed, for performance."""
        if self._n_refl[block_id] > 1:
            value, _, __ = self._smoother.multi_value_weight(
                self._normalised_x_values[block_id],
                self._normalised_y_values[block_id],
                self.value,
            )
        elif self._n_refl[block_id] == 1:
            value, _, __ = self._smoother.value_weight(
                self._normalised_x_values[block_id][0],
                self._normalised_y_values[block_id][0],
                self.value,
            )
            value = flex.double(1, value)
        else:
            value = flex.double([])
        return value


class SmoothScaleComponent3D(ScaleComponentBase, SmoothMixin):
    """Implementation of a 3D array-based smoothly varying scale factor.

    A 3d array of parameters is defined, and the scale factor at fractional
    coordinates is calculated as smoothly varying based on the distance to
    the nearby parameters as calculated in the GaussianSmoother3D. The
    initial values are passed as a 1D array, and shape is a 3-tuple
    indicating the number of parameters in each dimension."""

    null_parameter_value = 1.0

    def __init__(self, initial_values, shape, parameter_esds=None):
        assert len(initial_values) == (
            shape[0] * shape[1] * shape[2]
        ), """The
    shape information to initialise a 3D smoother is inconsistent with the
    length of the initial parameter list."""
        super(SmoothScaleComponent3D, self).__init__(initial_values, parameter_esds)
        self._n_x_params = shape[0]
        self._n_y_params = shape[1]
        self._n_z_params = shape[2]
        self._normalised_x_values = None
        self._normalised_y_values = None
        self._normalised_z_values = None

    def set_new_parameters(self, new_parameters, shape):
        """Set new parameters of a different length i.e. after batch handling"""
        assert len(new_parameters) == shape[0] * shape[1] * shape[2]
        self._parameters = new_parameters
        self._parameter_esds = None
        self._n_params = len(self._parameters)
        self._n_x_params = shape[0]
        self._n_y_params = shape[1]
        self._n_z_params = shape[2]

    @ScaleComponentBase.data.setter
    def data(self, data):
        assert set(data.keys()) == {"x", "y", "z"}, set(data.keys())
        self._data = data

    @property
    def n_x_params(self):
        """The number of parameters that parameterise the x-component."""
        return self._n_x_params

    @property
    def n_y_params(self):
        """The number of parameters that parameterise the y-component."""
        return self._n_y_params

    @property
    def n_z_params(self):
        """The number of parameters that parameterise the z-component."""
        return self._n_z_params

    @property
    def normalised_x_values(self):
        """The normalised coordinate values in the first dimension."""
        return self._normalised_x_values

    @property
    def normalised_y_values(self):
        """The normalised coordinate values in the second dimension."""
        return self._normalised_y_values

    @property
    def normalised_z_values(self):
        """The normalised coordinate values in the third dimension."""
        return self._normalised_z_values

    def update_reflection_data(self, selection=None, block_selections=None):
        """control access to setting all of reflection data at once"""
        self._normalised_x_values = []
        self._normalised_y_values = []
        self._normalised_z_values = []
        self._n_refl = []
        normalised_x_values = self.data["x"]
        normalised_y_values = self.data["y"]
        normalised_z_values = self.data["z"]
        if selection:
            normalised_x_values = normalised_x_values.select(selection)
            normalised_y_values = normalised_y_values.select(selection)
            normalised_z_values = normalised_z_values.select(selection)
        """Set the normalised coordinate values and configure the smoother."""
        normalised_x_values = normalised_x_values - flex.min(normalised_x_values)
        normalised_y_values = normalised_y_values - flex.min(normalised_y_values)
        normalised_z_values = normalised_z_values - flex.min(normalised_z_values)
        x_range = [
            floor(round(flex.min(normalised_x_values), 10)),
            ceil(round(flex.max(normalised_x_values), 10)),
        ]
        y_range = [
            floor(round(flex.min(normalised_y_values), 10)),
            ceil(round(flex.max(normalised_y_values), 10)),
        ]
        z_range = [
            floor(round(flex.min(normalised_z_values), 10)),
            ceil(round(flex.max(normalised_z_values), 10)),
        ]
        self._smoother = GaussianSmoother3D(
            x_range,
            self.nparam_to_val(self._n_x_params),
            y_range,
            self.nparam_to_val(self._n_y_params),
            z_range,
            self.nparam_to_val(self._n_z_params),
        )
        if block_selections:
            for i, sel in enumerate(block_selections):
                self._normalised_x_values.append(normalised_x_values.select(sel))
                self._normalised_y_values.append(normalised_y_values.select(sel))
                self._normalised_z_values.append(normalised_z_values.select(sel))
                self._n_refl.append(self._normalised_x_values[i].size())
        else:
            self._normalised_x_values.append(normalised_x_values)
            self._normalised_y_values.append(normalised_y_values)
            self._normalised_z_values.append(normalised_z_values)
            self._n_refl.append(normalised_x_values.size())

    def calculate_scales_and_derivatives(self, block_id=0):
        if self._n_refl[block_id] > 1:
            value, weight, sumweight = self._smoother.multi_value_weight(
                self._normalised_x_values[block_id],
                self._normalised_y_values[block_id],
                self._normalised_z_values[block_id],
                self.value,
            )
            inv_sw = 1.0 / sumweight
            dv_dp = row_multiply(weight, inv_sw)
        elif self._n_refl[block_id] == 1:
            value, weight, sumweight = self._smoother.value_weight(
                self._normalised_x_values[block_id][0],
                self._normalised_y_values[block_id][0],
                self._normalised_z_values[block_id][0],
                self.value,
            )
            dv_dp = sparse.matrix(1, weight.size)
            b = flex.double(weight.as_dense_vector() / sumweight)
            b.reshape(flex.grid(1, b.size()))
            dv_dp.assign_block(b, 0, 0)
            value = flex.double(1, value)
        else:
            return flex.double([]), sparse.matrix(0, 0)
        return value, dv_dp

    def calculate_scales(self, block_id=0):
        """"Only calculate the scales if needed, for performance."""
        if self._n_refl[block_id] > 1:
            value, _, __ = self._smoother.multi_value_weight(
                self._normalised_x_values[block_id],
                self._normalised_y_values[block_id],
                self._normalised_z_values[block_id],
                self.value,
            )
        elif self._n_refl[block_id] == 1:
            value, _, __ = self._smoother.value_weight(
                self._normalised_x_values[block_id][0],
                self._normalised_y_values[block_id][0],
                self._normalised_z_values[block_id][0],
                self.value,
            )
            value = flex.double(1, value)
        else:
            value = flex.double([])
        return value

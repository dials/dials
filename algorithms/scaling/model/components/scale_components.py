"""
Classes that define a component of a scaling model.

Each class holds the parameters and relevant data, as a list of
arrays, from which to calculate inverse scale factors and
derivatives.
The components are initialised without any data, which is added
by setting the data dict. In order to update the internal data
lists in order to calculate the scales and derivatives, the
update_reflection_data method should be called, which can optionally
be provided with selection arrays to split the data for blockwise/parallel
calculations.

The scaling algorithm makes use of the components in the following way.
First, the data for all 'suitable' reflections are added to the components.
Then, at different stages of the algorithm, selection lists are provided
to select a subset of this data (e.g. a small subset to prepare the
component for minimisation calculation, or a large subset for calculating
the scales for all reflections). The selection lists typically come from
the Ih_table datastructures so that the data in the components is split in
the same way as the data in the Ih_table datastructure.
"""
from __future__ import absolute_import, division, print_function
import abc
from dials.array_family import flex
from scitbx import sparse
from dials_scaling_ext import calculate_harmonic_tables_from_selections


class ScaleComponentBase(object):
    """
    Base scale component class.

    This defines an interface to access the parameters, the component
    of the inverse scale factor and it's derivatives with respect to
    the parameters. Scale components derived from the base class are
    designed to be instantiated by a ScalingModel class, by supplying
    an initial array of parameters and optionally the current estimated
    standard deviations. The relevant data from a reflection table is
    added later by a Scaler using the update_reflection_data method.
    This behaviour allows data to easily be added/changed after selecting
    subsets of the data.
    """

    __metaclass__ = abc.ABCMeta

    def __init__(self, initial_values, parameter_esds=None):
        """Set the initial parameter values, parameter esds and n_params."""
        self._parameters = initial_values
        self._parameter_esds = parameter_esds
        self._n_params = len(self._parameters)
        self._var_cov = None
        self._n_refl = []  # store as a list, to allow holding of data in blocks
        self._parameter_restraints = None
        self._data = {}

    @property
    def data(self):
        """
        Return a dictionary of reflection data relevant to the particular component.

        This is designed to be a dict of arrays which can be selected from when
        updating the component (i.e. selecting subsets).
        """
        return self._data

    @data.setter
    def data(self, data):
        self._data = data

    @property
    def parameter_restraints(self):
        """Restraint weights for the component parameters."""
        return self._parameter_restraints

    @parameter_restraints.setter
    def parameter_restraints(self, restraints):
        assert restraints.size() == self.parameters.size()
        self._parameter_restraints = restraints

    @property
    def n_params(self):
        """Get the number of parameters of the component (read-only)."""
        return self._n_params

    @property
    def parameters(self):
        """Parameters of the component."""
        return self._parameters

    @parameters.setter
    def parameters(self, new_parameters):
        assert len(new_parameters) == len(
            self._parameters
        ), """
attempting to set a new set of parameters of different length than previous
assignment: was %s, attempting %s""" % (
            len(self._parameters),
            len(new_parameters),
        )
        self._parameters = new_parameters

    @property
    def parameter_esds(self):
        """Return the estimated standard deviations of the parameters."""
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
        """Return the variance-covariance matrix of the parameters."""
        return self._var_cov

    @var_cov_matrix.setter
    def var_cov_matrix(self, var_cov):
        self._var_cov = var_cov

    @property
    def n_refl(self):
        """Return a list of the number of reflections in each block."""
        return self._n_refl

    @abc.abstractmethod
    def update_reflection_data(self, selection=None, block_selections=None):
        """
        Update the internal data arrays.

        Use the data stored in self.data, optionally with a selection array
        or list of selections, to populate a list of internal arrays e.g n_refl,
        normalised_values etc. to allow scale and derivative calculations. If no
        selection arrays are provided, the internal arrays will be lists
        containing one array/value, depending on the data type needed for
        derivative and scale calculation.

        Args:
            selection: A flex.bool selection array to select a subset of the
                internal data.
            block_selections (list): A list of flex.size_t arrays to select
                subsets of the internal data.
        """

    @abc.abstractmethod
    def calculate_scales_and_derivatives(self, block_id=0):
        """Calculate and return inverse scales and derivatives for a given block."""

    @abc.abstractmethod
    def calculate_scales(self, block_id=0):
        """Calculate and return inverse scales for a given block."""


class SingleScaleFactor(ScaleComponentBase):
    """
    A model component consisting of a single global scale parameter.

    The inverse scale factor for every reflection is the parameter
    value itself and the derivatives are therefore all 1.0.
    """

    null_parameter_value = 1.0

    def __init__(self, initial_values, parameter_esds=None):
        """Set the initial parameter values, parameter esds and n_params."""
        assert (
            len(initial_values) == 1
        ), """
This model component can only hold a single parameter."""
        super(SingleScaleFactor, self).__init__(initial_values, parameter_esds)

    @ScaleComponentBase.data.setter
    def data(self, data):
        """Set the data dict in the parent class."""
        assert set(data.keys()) == {"id"}, set(data.keys())
        self._data = data

    def update_reflection_data(self, selection=None, block_selections=None):
        """
        Update the internal n_refl list.

        Use the data stored in self.data, optionally with a boolean selection array
        or list of flex.size_t index selections, to make a list of n_refl (of length
        1 or len(block_selections)) by inspecting the size of the selection result,
        in order to allow scale and derivative calculations.

        Args:
            selection: Optional, a flex.bool selection array to select a subset of
                the internal data.
            block_selections (list): Optional, a list of flex.size_t arrays to
                select subsets of the internal data.
        """
        data = self.data["id"]
        if selection:
            self._n_refl = [data.select(selection).size()]
        elif block_selections:
            self._n_refl = [data.select(sel).size() for sel in block_selections]
        else:
            self._n_refl = [data.size()]

    def calculate_scales_and_derivatives(self, block_id=0):
        """Calculate and return inverse scales and derivatives for a given block."""
        scales = flex.double(self.n_refl[block_id], self._parameters[0])
        derivatives = sparse.matrix(self.n_refl[block_id], 1)
        for i in range(self.n_refl[block_id]):
            derivatives[i, 0] = 1.0
        return scales, derivatives

    def calculate_scales(self, block_id=0):
        """Calculate and return inverse scales for a given block."""
        return flex.double(self.n_refl[block_id], self._parameters[0])


class SingleBScaleFactor(ScaleComponentBase):
    """
    A model component for a single global B-factor parameter.

    The inverse scale factor for each reflection is given by
    S = exp(B/(2 * d^2)), the derivatives are S/(2 * d^2).
    """

    null_parameter_value = 0.0

    def __init__(self, initial_values, parameter_esds=None):
        """Set the initial parameter values, parameter esds and n_params."""
        super(SingleBScaleFactor, self).__init__(initial_values, parameter_esds)
        self._d_values = []

    @property
    def d_values(self):
        """Return a list of arrays of d-values associated with this component."""
        return self._d_values

    @ScaleComponentBase.data.setter
    def data(self, data):
        """Set the data dict in the parent class."""
        assert set(data.keys()) == {"id", "d"}, set(data.keys())
        self._data = data

    def update_reflection_data(self, selection=None, block_selections=None):
        """
        Update the internal n_refl and d_values lists.

        Use the data stored in self.data, optionally with a boolean selection array
        or list of flex.size_t index selections, to make a lists of n_refl and
        d_value arrays (of length 1 or len(block_selections)), in order to allow
        scale and derivative calculations.

        Args:
            selection: Optional, a flex.bool selection array to select a subset of
                the internal data.
            block_selections (list): Optional, a list of flex.size_t arrays to
                select subsets of the internal data.
        """
        data = self.data["d"]
        if selection:
            self._d_values = [data.select(selection)]
        elif block_selections:
            self._d_values = [data.select(sel) for sel in block_selections]
        else:
            self._d_values = [data]
        self._n_refl = [dvalues.size() for dvalues in self._d_values]

    def calculate_scales_and_derivatives(self, block_id=0):
        """Calculate and return inverse scales and derivatives for a given block."""
        d_squared = self._d_values[block_id] * self._d_values[block_id]
        scales = flex.exp(
            flex.double(self._n_refl[block_id], self._parameters[0]) / (2.0 * d_squared)
        )
        derivatives = sparse.matrix(self._n_refl[block_id], 1)
        for i in range(self._n_refl[block_id]):
            derivatives[i, 0] = scales[i] / (2.0 * d_squared[i])
        return scales, derivatives

    def calculate_scales(self, block_id=0):
        """Calculate and return inverse scales for a given block."""
        scales = flex.exp(
            flex.double(self._n_refl[block_id], self._parameters[0])
            / (2.0 * (self._d_values[block_id] * self._d_values[block_id]))
        )
        return scales


class SHScaleComponent(ScaleComponentBase):
    """
    A model component for a spherical harmonic absorption correction.

    This component uses a set of spherical harmonic functions to define
    an absorption surface for the crystal. A matrix of spherical harmonic
    coefficients for the data is stored in self._harmonic_values and is
    used to calculate the scales and derivatives.
    The scale is given by S = 1 + (sum_l sum_m Clm * Ylm) where Clm are
    the model parameters and Ylm are the spherical harmonic coefficients,
    the derivatives are then simply the coefficients Ylm.
    """

    null_parameter_value = 0.0
    coefficients_list = None  # shared class variable to reduce memory load

    def __init__(self, initial_values, parameter_esds=None):
        """Set the initial parameter values, parameter esds and n_params."""
        super(SHScaleComponent, self).__init__(initial_values, parameter_esds)
        self._harmonic_values = []
        self._matrices = []

    @property
    def harmonic_values(self):
        """Return the matrix of harmonic coefficients for the internal data."""
        return self._harmonic_values

    @property
    def sph_harm_table(self):
        """Return the matrix of the full harmonic coefficient for a reflection table."""
        return self._data["sph_harm_table"]

    @sph_harm_table.setter
    def sph_harm_table(self, sht):
        self._data["sph_harm_table"] = sht

    @ScaleComponentBase.data.setter
    def data(self, data):
        """Set the data dict in the parent class."""
        try:
            assert set(data.keys()) == {"s1_lookup", "s0_lookup"}, set(data.keys())
            self._mode = "memory"
        except AssertionError as e:
            assert set(data.keys()) == {"sph_harm_table"}, set(data.keys())
            self._mode = "speed"  # Note: only speedier for small datasets
        self._data = data

    def calculate_restraints(self):
        """Calculate residual and gradient restraints for the component."""
        residual = self.parameter_restraints * self._parameters * self._parameters
        gradient = 2.0 * self.parameter_restraints * self._parameters
        return residual, gradient

    def calculate_jacobian_restraints(self):
        """Calculate residual and jacobian restraints for the component."""
        jacobian = sparse.matrix(self.n_params, self.n_params)
        for i in range(self.n_params):
            jacobian[i, i] = 1.0
        return self._parameters, jacobian, self._parameter_restraints

    def update_reflection_data(self, selection=None, block_selections=None):
        """
        Update the internal n_refl and harmonic_values lists.

        Use the harmonic values matrix stored in self.data, optionally with a
        boolean selection array or list of flex.size_t index selections, to make
        lists of n_refl and harmonic_value arrays (of length 1 or
        len(block_selections)), in order to allow scale and derivative calculations.

        Args:
            selection: Optional, a flex.bool selection array to select a subset of
                the internal data.
            block_selections (list): Optional, a list of flex.size_t arrays to
                select subsets of the internal data.
        """
        if self._mode == "speed":
            self._update_reflection_data_speedmode(selection, block_selections)
        elif self._mode == "memory":
            self._update_reflection_data_memorymode(selection, block_selections)
        else:
            raise ValueError

    def _update_reflection_data_memorymode(self, selection=None, block_selections=None):
        if len(self.coefficients_list) != self.n_params:
            self.coefficients_list = self.coefficients_list[0 : self.n_params]
            # modify only for this instance, only needs to be done once per instance.
        if selection:
            n0 = self.data["s0_lookup"].select(selection)
            n1 = self.data["s1_lookup"].select(selection)
            values, matrix = calculate_harmonic_tables_from_selections(
                n0, n1, self.coefficients_list
            )
            self._harmonic_values = [values]
            self._matrices = [matrix]
        elif block_selections:
            self._harmonic_values = []
            self._matrices = []
            for sel in block_selections:
                n0 = self.data["s0_lookup"].select(sel)
                n1 = self.data["s1_lookup"].select(sel)
                values, matrix = calculate_harmonic_tables_from_selections(
                    n0, n1, self.coefficients_list
                )
                self._harmonic_values.append(values)
                self._matrices.append(matrix)
        else:
            n0 = self.data["s0_lookup"]
            n1 = self.data["s1_lookup"]
            values, matrix = calculate_harmonic_tables_from_selections(
                n0, n1, self.coefficients_list
            )
            self._harmonic_values = [values]
            self._matrices = [matrix]
        self._n_refl = [val[0].size() for val in self._harmonic_values]

    def _update_reflection_data_speedmode(self, selection=None, block_selections=None):
        if selection:
            sel_sph_harm_table = self.data["sph_harm_table"].select_columns(
                selection.iselection()
            )
            self._harmonic_values = [sel_sph_harm_table.transpose()]
        elif block_selections:
            self._harmonic_values = []
            for sel in block_selections:
                block_sph_harm_table = self.data["sph_harm_table"].select_columns(sel)
                self._harmonic_values.append(block_sph_harm_table.transpose())
        else:
            self._harmonic_values = [self.data["sph_harm_table"].transpose()]
        self._n_refl = [val.n_rows for val in self._harmonic_values]

    def calculate_scales(self, block_id=0):
        """Calculate and return inverse scales for a given block."""
        if self._mode == "speed":
            return self._calculate_scales_and_derivatives_speedmode(
                block_id, derivatives=False
            )
        elif self._mode == "memory":
            return self._calculate_scales_and_derivatives_memorymode(
                block_id, derivatives=False
            )

    def calculate_scales_and_derivatives(self, block_id=0):
        """Calculate and return inverse scales and derivatives for a given block."""
        if self._mode == "speed":
            return self._calculate_scales_and_derivatives_speedmode(block_id)
        elif self._mode == "memory":
            return self._calculate_scales_and_derivatives_memorymode(block_id)

    def _calculate_scales_and_derivatives_speedmode(self, block_id, derivatives=True):
        abs_scale = flex.double(
            self._harmonic_values[block_id].n_rows, 1.0
        )  # Unity term
        for i, col in enumerate(self._harmonic_values[block_id].cols()):
            abs_scale += flex.double(col.as_dense_vector() * self._parameters[i])
        if derivatives:
            return abs_scale, self._harmonic_values[block_id]
        return abs_scale

    def _calculate_scales_and_derivatives_memorymode(self, block_id, derivatives=True):
        abs_scale = flex.double(
            self._harmonic_values[block_id][0].size(), 1.0
        )  # Unity term
        for i, arr in enumerate(
            self._harmonic_values[block_id]
        ):  # iterate over a list of arrays
            abs_scale += arr * self._parameters[i]
        if derivatives:
            return abs_scale, self._matrices[block_id]
        return abs_scale

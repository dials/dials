"""
Classes to initialise a 'parameter manager', to indicate to a
refiner which components of the model are to be refined.
"""
from __future__ import absolute_import, division, print_function

import logging
from collections import OrderedDict

from dials.array_family import flex
import six

logger = logging.getLogger("dials")


class active_parameter_manager(object):
    """
    Class to manage the current active parameters during minimisation.

    The parameter manager is initialised with components - a dict of ('name' : obj)
    pairs, and a selection list, which is a list of 'name' values to be minimised
    in this cycle. obj is the component of a model, and this code requires that
    obj has the attribute 'parameters'.
    This class stores a reference to the components to be minimised, alongside
    some bookkeeping to allow selection of the parameters for an individual
    component.
    """

    def __init__(self, components, selection_list):
        self.x = flex.double([])
        self.components = OrderedDict()
        self.derivatives = None
        self.var_cov_matrix = None
        self.components_list = []  # just a list of the component names
        n_cumul_params = 0
        for component, obj in six.iteritems(components):
            if component in selection_list:
                assert hasattr(
                    obj, "parameters"
                ), """component object must have the
          attribute 'parameters' for access to the component parameters."""
                self.x.extend(obj.free_parameters)
                n_params = len(obj.free_parameters)
                self.components.update(
                    {
                        component: {
                            "object": obj,
                            "n_params": n_params,
                            "start_idx": n_cumul_params,
                            "end_idx": len(self.x),
                        }
                    }
                )
                n_cumul_params += n_params
        self.n_active_params = len(self.x)
        for comp in self.components:
            self.components_list.extend([comp])

    def select_parameters(self, component):
        """Select the subset of self.x corresponding to the component (a string)."""
        start_idx = self.components[component]["start_idx"]
        end_idx = self.components[component]["end_idx"]
        return self.x[start_idx:end_idx]

    def set_param_vals(self, x):
        """Set method for refinement engine access."""
        self.x = x
        for component in self.components:
            component_obj = self.components[component]["object"]
            component_obj.free_parameters = self.select_parameters(component)

    def get_param_vals(self):
        """Get method for refinement engine access."""
        return self.x

    def calculate_model_state_uncertainties(self, var_cov):
        """Set var_cov matrices for each component to allow later calculation
        of errors."""
        i = 0
        for component in self.components.values():
            n = component["n_params"]
            sub = var_cov.matrix_copy_block(i, i, n, n)
            component["object"].var_cov_matrix = sub
            i += n
        self.var_cov_matrix = var_cov

    def set_param_esds(self, esds):
        """Set the estimated standard deviations of the model components."""
        for component in self.components.values():
            start_idx = component["start_idx"]
            end_idx = component["end_idx"]
            component["object"].free_parameter_esds = esds[start_idx:end_idx]


class multi_active_parameter_manager(object):
    """
    Parameter manager to manage the current active parameters during minimisation
    for multiple datasets that are being minimised simultaneously.

    Initialise with two lists of components and selections, each item of which
    is used to initialise an active parameter manager of type apm_class.
    """

    def __init__(self, components_list, selection_lists, apm_class):
        self.x = flex.double([])
        self.derivatives = None
        self.components_list = []  # A list of the component names.
        self.apm_list = []
        self.apm_data = OrderedDict()
        all_same_components = False
        if all(i == selection_lists[0] for i in selection_lists[1:]):
            logger.info(
                "Components to be refined in this cycle for all datasets: %s",
                "".join(str(i) + ", " for i in selection_lists[0]).rstrip(", "),
            )
            all_same_components = True
        for j, (components, selection_list) in enumerate(
            zip(components_list, selection_lists)
        ):
            self.apm_list.append(apm_class(components, selection_list))
            if not all_same_components:
                logger.info(
                    "Components to be refined in this cycle for datasest %s: %s",
                    j,
                    ",".join(i for i in self.apm_list[j].components_list),
                )
        n_cumul_params = 0
        for i, apm in enumerate(self.apm_list):
            self.x.extend(apm.x)
            n_params = apm.n_active_params
            self.apm_data.update(
                {i: {"start_idx": n_cumul_params, "end_idx": len(self.x)}}
            )
            n_cumul_params += n_params
        self.n_active_params = len(self.x)

        for apm in self.apm_list:
            for comp in apm.components:
                self.components_list.extend([comp])

    def select_parameters(self, apm_number):
        """Select the subset of self.x corresponding to the apm number."""
        apm_data = self.apm_data[apm_number]
        return self.x[apm_data["start_idx"] : apm_data["end_idx"]]

    def set_param_vals(self, x):
        """Set method for refinement engine access."""
        self.x = x
        for i, single_apm in enumerate(self.apm_list):
            single_apm.set_param_vals(self.select_parameters(i))

    def get_param_vals(self):
        """Get method for refinement engine access."""
        return self.x

    def set_param_esds(self, esds):
        """Set the estimated standard deviations of the parameters."""
        for i, apm in enumerate(self.apm_list):
            apm_data = self.apm_data[i]
            apm.set_param_esds(esds[apm_data["start_idx"] : apm_data["end_idx"]])

    def calculate_model_state_uncertainties(self, var_cov):
        """Set var_cov matrices for each component, to allow later calculation
        of errors."""
        i = 0
        for apm in self.apm_list:
            n = apm.n_active_params
            sub = var_cov.matrix_copy_block(i, i, n, n)
            apm.calculate_model_state_uncertainties(sub)
            i += n


class ParameterManagerGenerator(object):
    """
    Class to generate multi-dataset parameter managers for minimisation.

    Handles the case of concurrent parameter minimisation (one parameter
    manager generated) and the case of consecutive parameter minimisation
    (several parameter managers generated, one for each minimisation, depending
    on the data_manager.consecutive_refinement_order property).
    """

    def __init__(self, data_managers, apm_type, mode="concurrent"):
        if mode not in ["concurrent", "consecutive"]:
            raise ValueError(
                "Bad value for refinement order mode: %s, expected %s"
                % (mode, " or ".join(["concurrent", "consecutive"]))
            )
        self.data_managers = data_managers
        self.apm_type = apm_type
        self.mode = mode
        self.param_lists = [None] * len(data_managers)
        if self.mode == "concurrent":
            for i, data_manager in enumerate(self.data_managers):
                self.param_lists[i] = [param for param in data_manager.components]
        else:  # mode=consecutive
            # Generate nested list indicating the names of active parameters
            # e.g consecutive_order for class is [["a", "b"], ["c"]],
            # data_man has components ['a'], return [["a"]]
            # data_man has components ['a', 'c'], return [["a"], ["c"]]
            for i, data_manager in enumerate(self.data_managers):
                ind_param_list = []
                for cycle in data_manager.consecutive_refinement_order:
                    corrlist = [
                        corr for corr in cycle if corr in data_manager.components
                    ]
                    if corrlist:
                        ind_param_list.append(corrlist)
                self.param_lists[i] = ind_param_list

    def parameter_managers(self):
        """Generate the parameter managers for each cycle of refinement."""
        if self.mode == "concurrent":
            return self._parameter_managers_concurrent()
        return self._parameter_managers_consecutive()

    def _parameter_managers_concurrent(self):
        components = [s.components for s in self.data_managers]
        return [
            multi_active_parameter_manager(components, self.param_lists, self.apm_type)
        ]

    def _parameter_managers_consecutive(self):
        components = [s.components for s in self.data_managers]
        while any(self.param_lists[i] for i, _ in enumerate(self.data_managers)):
            params = []
            for param_list in self.param_lists:
                if param_list:
                    params.append(param_list.pop(0))
                else:
                    params.append([])
            yield multi_active_parameter_manager(components, params, self.apm_type)

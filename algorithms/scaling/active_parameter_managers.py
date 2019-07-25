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
        self.curvatures = None
        self.var_cov_matrix = None
        self.components_list = []  # just a list of the component names
        n_cumul_params = 0
        for component, obj in six.iteritems(components):
            if component in selection_list:
                assert hasattr(
                    obj, "parameters"
                ), """component object must have the
          attribute 'parameters' for access to the component parameters."""
                self.x.extend(obj.parameters)
                n_params = len(obj.parameters)
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
        # logger.info('Components to be refined in this cycle: %s \n',
        #  ''.join(str(i)+', ' for i in self.components_list).rstrip(', '))

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
            component_obj.parameters = self.select_parameters(component)

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
            component["object"].parameter_esds = esds[start_idx:end_idx]


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
        self.curvatures = None
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
                    "".join(str(i) + ", " for i in components).rstrip(", "),
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
        logger.info(
            "\nConfigured a parameter manager for %s datasets.\n", len(self.apm_list)
        )

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


class ConcurrentAPMFactory(object):
    """
    Factory to correctly set up a single/multi active parameter manager for
    concurrent scaling of all model components of a general data_manager
    (or multiple data managers). This extracts the name of the components from
    each data_manager and creates a list to pass on to the apm_type specified
    (e.g a subclass of active_parameter_manager). mode=single/multi returns a
    single/multi active parameter manager.

    Data managers is a list of objects which have a components attribute.
    """

    def __init__(self, data_managers, apm_type, multi_mode=True):
        # One can optionally set multi_mode=False to force a single_apm to be created
        # for only one dataset - however all scaling methods are designed to iterate
        # over a multi_apm with an apm_list property.
        self.data_managers = data_managers
        self.apm = None
        self.multi_mode = multi_mode
        if len(data_managers) > 1:
            self.multi_mode = True
        self.param_lists = []
        self.create_active_list(apm_type)
        self.n_cycles = 1

    def create_active_list(self, apm_type):
        """Return a list indicating the names of active parameters."""

        if not self.multi_mode:
            param_name = []
            for param in self.data_managers[0].components:
                param_name.append(str(param))
            if not param_name:
                raise ValueError(
                    "No model components have been chosen, aborting process."
                )
            self.param_lists = param_name
            self.apm = apm_type(self.data_managers[0].components, self.param_lists)

        else:
            for data_manager in self.data_managers:
                param_name = []
                for param in data_manager.components:
                    param_name.append(str(param))
                if not param_name:
                    raise ValueError(
                        "No model components have been chosen, aborting process."
                    )
                self.param_lists.append(param_name)
            components = [i.components for i in self.data_managers]
            self.apm = multi_active_parameter_manager(
                components, self.param_lists, apm_type
            )

    def make_next_apm(self):
        """Method to call to return the apm."""
        return self.apm


class ConsecutiveAPMFactory(object):
    """
    Factory to correctly set up a nested list structure to pass to
    single/multi active parameter managers for consecutive scaling of the
    model components of a general data_manager (or multiple data managers).

    Upon calling make_next_apm, the first list element for each dataset is used to
    initialise the apm_type specified (e.g a subclass of active_parameter_manager)
    and then removed from the list structure.
    make_next_apm can be called n_cycles times.
    mode=single/multi returns a single/mutli active parameter manager.
    """

    def __init__(self, data_managers, apm_type, multi_mode=True):
        # One can optionally set multi_mode=False to force a single_apm to be created
        # for only one dataset - however all scaling methods are designed to iterate
        # over a multi_apm with an apm_list property.
        self.data_managers = data_managers
        self.multi_mode = multi_mode
        if len(data_managers) > 1:
            self.multi_mode = True
        self.apm_type = apm_type
        self.param_lists = []
        self.n_cycles = None
        self.create_consecutive_list()

    def create_consecutive_list(self):
        """Return a list indicating the names of active parameters."""

        if not self.multi_mode:
            for cycle in self.data_managers[0].consecutive_refinement_order:
                corrlist = []
                for corr in cycle:
                    if corr in self.data_managers[0].components:
                        corrlist.append(corr)
                self.param_lists.append(corrlist)
            self.n_cycles = sum([1 for i in self.param_lists if i])

        else:
            for data_manager in self.data_managers:
                ind_param_list = []
                for cycle in data_manager.consecutive_refinement_order:
                    corrlist = []
                    for corr in cycle:
                        if corr in data_manager.components:
                            corrlist.append(corr)
                    ind_param_list.append(corrlist)
                self.param_lists.append(ind_param_list)
            # now need to calculate the max number of cycles needed across all data_managers
            is_cycle_active = []
            for p_list in self.param_lists:
                for i, cycle in enumerate(p_list):
                    if cycle:
                        is_cycle_active.append(i)
            self.n_cycles = len(set(is_cycle_active))
            # now make sure all lists are same length
            max_len = max([len(i) for i in self.param_lists])
            for p_list in self.param_lists:
                for _ in range(max_len - len(p_list)):
                    p_list.append([])

    def make_next_apm(self):
        """Generate a valid apm for minimisation (contains some active parameters,
        but not necessarily for all datasets)."""

        if not self.multi_mode:
            apm = self.apm_type(self.data_managers[0].components, self.param_lists[0])
            self.param_lists = self.param_lists[1:]
        else:
            params = []
            for i in range(0, len(self.param_lists)):
                params.append(self.param_lists[i][0])
                self.param_lists[i] = self.param_lists[i][1:]  # remove first element
            components = [i.components for i in self.data_managers]
            apm = multi_active_parameter_manager(components, params, self.apm_type)
        if not apm.components_list:  # no active parameters, so iterate
            apm = self.make_next_apm()
        return apm

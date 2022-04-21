"""
Classes to initialise a 'parameter manager', to indicate to a
refiner which components of the model are to be refined.
"""

from __future__ import annotations

import logging

from scitbx import sparse

from dials.array_family import flex

logger = logging.getLogger("dials")


class active_parameter_manager:
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
        # self.target = target
        self.x = flex.double([])
        self.components = {}
        self.derivatives = None
        self.var_cov_matrix = None
        self.components_list = []  # just a list of the component names
        n_cumul_params = 0
        for component, obj in components.items():
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

    def get_param_names(self):
        l = []
        for name, data in self.components.items():
            l.extend([name + f"_param{j}" for j in range(data["n_params"])])
        return l

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


class TargetInterface:
    def compute_functional_gradients(self, block):
        return self.target.compute_functional_gradients(block)

    def compute_restraints_functional_gradients(self, block):
        return self.target.compute_restraints_functional_gradients(block)

    def compute_residuals_and_gradients(self, block):
        return self.target.compute_residuals_and_gradients(block)

    def compute_restraints_residuals_and_gradients(self, block):
        return self.target.compute_restraints_residuals_and_gradients(block)

    def compute_residuals(self, block):
        return self.target.compute_residuals(block)


class multi_active_parameter_manager(TargetInterface):
    """
    Parameter manager to manage the current active parameters during minimisation
    for multiple datasets that are being minimised simultaneously.

    Initialise with two lists of components and selections, each item of which
    is used to initialise an active parameter manager of type apm_class.
    """

    def __init__(self, target, components_list, selection_lists, apm_class):
        self.target = target
        self.x = flex.double([])
        self.derivatives = None
        self.components_list = []  # A list of the component names.
        self.apm_list = []
        self.apm_data = {}
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
                    "Components to be refined in this cycle for dataset %s: %s",
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

    def get_param_names(self):
        l = []
        for i, apm in enumerate(self.apm_list):
            for name, data in apm.components.items():
                l.extend([name + f"_expt{i}_param{j}" for j in range(data["n_params"])])
        return l

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


class shared_active_parameter_manager(multi_active_parameter_manager):

    """Class to enforce sharing of model components.

    Intercept calls to a multi_apm, to override set_params calls and manage
    reshaping of gradient/jacobian during minimisation.
    """

    def __init__(self, target, components_list, selection_lists, apm_class, shared):
        super().__init__(target, components_list, selection_lists, apm_class)
        n_unique_params = 0
        found_initial_shared = False
        # first loop over to work out how many unique parameters overall
        for i, apm in enumerate(self.apm_list):
            for name, comp in apm.components.items():
                if name != shared:
                    n_unique_params += comp["n_params"]
                elif not found_initial_shared:
                    n_unique_params += comp["n_params"]
                    found_initial_shared = True

        self.reducing_matrix = sparse.matrix(self.n_active_params, n_unique_params)
        # now loop through to work out reducing matrix
        # also need to update apm_data with a size_t selection
        found_initial_shared = False
        shared_params = (0, 0)
        cumul_params = 0
        unique_parameters = []
        for i, apm in enumerate(self.apm_list):
            apm_sel = flex.size_t()
            # ^ needs to be size_t to make sure selected in correct order
            for name, comp in apm.components.items():
                start_col_idx = cumul_params
                indiv_start = self.apm_data[i]["start_idx"] + comp["start_idx"]
                indiv_end = self.apm_data[i]["start_idx"] + comp["end_idx"]
                if name != shared:
                    apm_sel.extend(
                        flex.size_t(
                            range(cumul_params, cumul_params + comp["n_params"])
                        )
                    )
                    for n, j in enumerate(range(indiv_start, indiv_end)):
                        self.reducing_matrix[j, start_col_idx + n] = 1
                        unique_parameters.append(j)
                    cumul_params += comp["n_params"]  #
                elif not found_initial_shared:
                    apm_sel.extend(
                        flex.size_t(
                            range(cumul_params, cumul_params + comp["n_params"])
                        )
                    )
                    for n, j in enumerate(range(indiv_start, indiv_end)):
                        self.reducing_matrix[j, start_col_idx + n] = 1
                        unique_parameters.append(j)
                    shared_start_column = start_col_idx
                    found_initial_shared = True
                    shared_params = (indiv_start, indiv_end)
                    cumul_params += comp["n_params"]
                else:  # name == shared and n_shared_found > 0:
                    assert (
                        shared_params[1] - shared_params[0] == indiv_end - indiv_start
                    )
                    apm_sel.extend(
                        flex.size_t(range(shared_params[0], shared_params[1]))
                    )
                    for n, j in enumerate(range(indiv_start, indiv_end)):
                        self.reducing_matrix[j, shared_start_column + n] = 1
            self.apm_data[i]["apm_sel"] = apm_sel

        # also need to reduce self.x to the size of the new params
        self.x = self.x.select(flex.size_t(unique_parameters))

    def select_parameters(self, apm_number):
        """Select the subset of self.x corresponding to the apm number."""
        apm_data = self.apm_data[apm_number]
        sel = apm_data["apm_sel"]
        return self.x.select(sel)

    def compute_functional_gradients(self, block):
        f, g = self.target.compute_functional_gradients(block)
        return f, g * self.reducing_matrix

    def compute_restraints_functional_gradients(self, block):
        res = self.target.compute_restraints_functional_gradients(block)
        if res is not None:
            return res[0], res[1] * self.reducing_matrix
        return res

    def compute_residuals_and_gradients(self, block):
        r, j, w = self.target.compute_residuals_and_gradients(block)
        return r, j * self.reducing_matrix, w

    def compute_restraints_residuals_and_gradients(self, block):
        res = self.target.compute_restraints_residuals_and_gradients(block)
        if res is not None:
            return res[0], res[1] * self.reducing_matrix, res[2]
        return res

    def set_param_esds(self, esds):
        """Set the estimated standard deviations of the parameters."""
        for apm, apm_data in zip(self.apm_list, self.apm_data.values()):
            apm.set_param_esds(esds.select(apm_data["apm_sel"]))

    def calculate_model_state_uncertainties(self, var_cov):
        """Set var_cov matrices for each component, to allow later calculation
        of errors."""
        for i, apm in enumerate(self.apm_list):
            sub_var_cov = sparse.matrix(apm.n_active_params, apm.n_active_params)
            n_this = 0
            for comp in apm.components.values():
                n = comp["n_params"]
                start_idx = self.apm_data[i]["apm_sel"][n_this]
                sub = var_cov.matrix_copy_block(start_idx, start_idx, n, n)
                sub_var_cov.assign_block(sub, n_this, n_this)
                n_this += n
            apm.calculate_model_state_uncertainties(sub_var_cov.as_dense_matrix())


class ParameterManagerGenerator:
    """
    Class to generate multi-dataset parameter managers for minimisation.

    Handles the case of concurrent parameter minimisation (one parameter
    manager generated) and the case of consecutive parameter minimisation
    (several parameter managers generated, one for each minimisation, depending
    on the data_manager.consecutive_refinement_order property).
    """

    def __init__(self, data_managers, apm_type, target, mode="concurrent", shared=None):
        if mode not in ["concurrent", "consecutive"]:
            raise ValueError(
                "Bad value for refinement order mode: %s, expected %s"
                % (mode, " or ".join(["concurrent", "consecutive"]))
            )
        self.target = target
        self.data_managers = data_managers
        self.apm_type = apm_type
        self.mode = mode
        self.shared = shared
        self.param_lists = [None] * len(data_managers)
        if self.mode == "concurrent":
            for i, data_manager in enumerate(self.data_managers):
                components = list(data_manager.components)
                for f in data_manager.fixed_components:
                    try:
                        components.remove(f)
                    except ValueError:
                        continue
                self.param_lists[i] = components
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
                    for f in data_manager.fixed_components:
                        try:
                            corrlist.remove(f)
                        except ValueError:
                            continue
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
        if self.shared:
            return [
                shared_active_parameter_manager(
                    self.target,
                    components,
                    self.param_lists,
                    self.apm_type,
                    self.shared,
                )
            ]
        return [
            multi_active_parameter_manager(
                self.target, components, self.param_lists, self.apm_type
            )
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
            if self.shared:
                yield shared_active_parameter_manager(
                    self.target, components, params, self.apm_type, self.shared
                )
            else:
                yield multi_active_parameter_manager(
                    self.target, components, params, self.apm_type
                )

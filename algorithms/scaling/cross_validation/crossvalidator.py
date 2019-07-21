"""
This module defines an abstract CrossValidator and an implementation of a
cross validator for dials.scale
"""
from __future__ import absolute_import, division, print_function

import abc
import itertools
from copy import deepcopy

from dials.algorithms.scaling.observers import register_merging_stats_observers
from libtbx.table_utils import simple_table
from scitbx.array_family import flex
import six


class CrossValidator(object):

    """Abstract class defining common methods for cross validation and methods
    that must be implemented for concrete implementations"""

    __metaclass__ = abc.ABCMeta

    # metadata needed when constructing the results table
    results_metadata = {
        "names": [],  # the names of the 'results' that should be printed in the
        # results table upon completion
        "indices_to_monitor": [],  # these indices of the above list will be
        # monitored to see which config gives the best values
        "best_criterion": [],  # define 'best' for the above indices ('min' or 'max')
    }

    def __init__(self, experiments, reflections):
        self.experiments = experiments
        self.reflections = reflections
        self.results_dict = {}

    @abc.abstractmethod
    def run_script(self, params, config_no):
        """Run the appropriate command line script with the params, get the
        free/work set results and add to the results dict. Indicate the
        configuration number being run."""

    @abc.abstractmethod
    def get_results_from_script(self, script):
        """Return the work/free results list from the command line script object"""

    @abc.abstractmethod
    def get_parameter_type(self, name):
        """Find the parameter type for a discreet phil option - bool or choice."""

    @abc.abstractmethod
    def set_parameter(self, params, name, val):
        """Find the name in the params scope extract and set it to the val"""

    @abc.abstractmethod
    def set_free_set_offset(self, params, n):
        """Set the free set offset in the correct place in the scope"""

    @abc.abstractmethod
    def get_free_set_percentage(self, params):
        """Inspect the free set percentage in the correct place in the scope"""

    @staticmethod
    def _avg_sd_from_list(lst):
        """simple function to get average and standard deviation"""
        arr = flex.double(lst)
        avg = round(flex.mean(arr), 5)
        std = round(arr.standard_deviation_of_the_sample(), 5)
        return avg, std

    def create_results_dict(self, n_options):
        """Create a results dict of the correct size for filling later"""
        if n_options == 1:
            self.results_dict[0] = {"configuration": ["user"]}
            for name in self.results_metadata["names"]:
                self.results_dict[0][name] = []
        else:
            self.results_dict = {}
            for i in range(n_options):
                self.results_dict[i] = {"configuration": []}
                for name in self.results_metadata["names"]:
                    self.results_dict[i][name] = []

    def set_results_dict_configuration(self, keys, values):
        """Add configuration information to the results dict"""
        assert len(keys) == len(values)
        for i, v in enumerate(itertools.product(*values)):
            e = dict(zip(keys, v))
            for k, val in six.iteritems(e):
                self.results_dict[i]["configuration"].append(str(k) + "=" + str(val))

    def add_results_to_results_dict(self, config_no, results):
        """Add the results to the correct place in the dict"""
        assert len(results) == len(self.results_metadata["names"])
        for name, result in zip(self.results_metadata["names"], results):
            self.results_dict[config_no][name].append(result)

    def interpret_results(self):
        """Inspect the data in results_dict, make a nice table with the mean and
        average over many attempts and indicate the 'best' option"""
        rows = []
        headers = ["option", ""] + self.results_metadata["names"]
        monitored_values = []

        # Construct the rows, using the metadata from the results dict
        for v in self.results_dict.values():
            config_str = " ".join(v["configuration"])
            vals, stds = [], []
            for i, name in enumerate(self.results_metadata["names"]):
                val, std = self._avg_sd_from_list(v[name])
                vals.append(val)
                stds.append(std)
                if i in self.results_metadata["indices_to_monitor"]:
                    monitored_values.append(val)
            rows.append([config_str, "mean"] + [str(i) for i in vals])
            rows.append(["", "std dev"] + [str(i) for i in stds])

        # Now go through monitored values, finding the best and adding a '*'
        n_monitored = len(self.results_metadata["indices_to_monitor"])
        for i in range(n_monitored):
            vals = monitored_values[i::n_monitored]
            if self.results_metadata["best_criterion"][i] == "max":
                best_idx = vals.index(max(vals)) * 2  # *2 to skip std rows
            elif self.results_metadata["best_criterion"][i] == "min":
                best_idx = vals.index(min(vals)) * 2  # *2 to skip std rows
            rows[best_idx][self.results_metadata["indices_to_monitor"][i] + 2] += "*"
            # line above, 2 is to offset first two columns in table

        return simple_table(rows, headers)


class DialsScaleCrossValidator(CrossValidator):

    """An implementation of the CrossValidator for running dials.scale"""

    results_metadata = {  # metadata used when constructing the results table
        "names": [
            "work Rpim",
            "free Rpim",
            "Rpim gap",
            "work CC1/2",
            "free CC1/2",
            "CC1/2 gap",
            "work w.a.-CC1/2",
            "free w.a.-CC1/2",
            "w.avg.-CC1/2 gap",
        ],
        "indices_to_monitor": [
            1,
            2,
            4,
            5,
            7,
            8,
        ],  # these indices of the above list will be
        # monitored to see which config gives the best values
        "best_criterion": ["min", "min", "max", "min", "max", "min"],
    }

    def get_results_from_script(self, script):
        """Return the work/free results list from the command line script object"""
        result = script.scaler.final_rmsds
        return result

    def get_parameter_type(self, name):
        """Find the parameter type for a discreet phil option - bool or choice."""
        # Note - ideally could inspect the phil_scope?
        if name in ["outlier_rejection", "model"]:
            return "choice"
        elif name in [
            "absorption_term",
            "decay_term",
            "scale_term",
            "modulation_term",
            "optimise_errors",
            "full_matrix",
            "concurrent",
            "target_cycle",
        ]:
            return "bool"
        elif name in [
            "lmax",
            "n_modulation_bins",
            "n_resolution_bins",
            "n_absorption_bins",
        ]:
            return "int"
        elif name in [
            "surface_weight",
            "scale_interval",
            "decay_interval",
            "d_min",
            "d_max",
            "outlier_zmax",
        ]:
            return "float"

    def set_parameter(self, params, name, val):
        """Find the name in the params scope extract and set it to the val"""
        # Note: must be a better way to do this?
        if name in [
            "lmax",
            "n_modulation_bins",
            "n_resolution_bins",
            "n_absorption_bins",
        ]:
            params.parameterisation.__setattr__(name, int(val))  # convert float to int
        elif name in [
            "scale_term",
            "scale_interval",
            "decay_term",
            "decay_interval",
            "absorption_term",
            "surface_weight",
            "modulation_term",
        ]:
            params.parameterisation.__setattr__(name, val)
        elif name in ["optimise_errors"]:
            params.weighting.__setattr__(name, val)
        elif name in ["d_min", "d_max"]:  # But what about biasing by n_refl?
            params.cut_data.__setattr__(name, val)
        elif name in [
            "target_cycle",
            "concurrent",
            "full_matrix",
            "outlier_zmax",
            "outlier_rejection",
        ]:
            params.scaling_options.__setattr__(name, val)
        elif name in ["model"]:
            params.__setattr__(name, val)
        else:
            assert 0, "Unable to set chosen attribute " + str(name) + "=" + str(val)
        return params

    def set_free_set_offset(self, params, n):
        """Set the free set offset in the correct place in the scope"""
        params.scaling_options.free_set_offset = n
        return params

    def get_free_set_percentage(self, params):
        """Inspect the free set percentage in the correct place in the scope"""
        return params.scaling_options.free_set_percentage

    def run_script(self, params, config_no):
        """Run the scaling script with the params, get the free/work set results
        and add to the results dict"""
        from dials.command_line.scale import Script

        params.scaling_options.__setattr__("use_free_set", True)
        script = Script(
            params,
            experiments=deepcopy(self.experiments),
            reflections=deepcopy(self.reflections),
        )
        register_merging_stats_observers(script)
        script.run()
        results = self.get_results_from_script(script)
        self.add_results_to_results_dict(config_no, results)

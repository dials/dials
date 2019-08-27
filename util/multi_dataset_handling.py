# coding: utf-8
"""
Module of functions for handling operations on lists of reflection tables
and experiment lists.
"""

from __future__ import absolute_import, division, print_function

import logging
import uuid

logger = logging.getLogger("dials")

import iotbx.phil

phil_scope = iotbx.phil.parse(
    """
  dataset_selection {
    use_datasets = None
      .type = strings
      .help = "Choose a subset of datasets, based on the dataset id (as defined
               in the reflection table), to use from a multi-dataset input."
      .expert_level = 2
    exclude_datasets = None
      .type = strings
      .help = "Choose a subset of datasets, based on the dataset id (as defined
               in the reflection table), to exclude from a multi-dataset input."
      .expert_level = 2
  }
"""
)


def generate_experiment_identifiers(experiments, identifier_type="uuid"):
    """Generate unique identifiers for each experiment."""
    if identifier_type == "uuid":
        for expt in experiments:
            expt.identifier = str(uuid.uuid4())
    elif identifier_type == "timestamp":
        pass


def split_reflection_tables_on_ids(reflection_tables):
    """"Split a list of multi-dataset reflection tables, selecting on id.

    Note that reflection table entries with -1 are removed, as they are not
    assigned to any experiment.

    Args:
        reflection_tables (list): A list of reflection tables (each of which may include
            multiple datasets).

    Returns:
        (list): A list of single-dataset reflection tables.
    """
    single_reflection_tables = []
    for table in reflection_tables:
        if len(set(table["id"]).difference({-1})) > 1:
            ##FIXME fix split_by_experiment_id so that don't need to filter
            # unindexed reflections here to get rid of id = -1
            if -1 in table["id"]:
                table = table.select(table["id"] != -1)
            #  split on id and preserve experiment_identifiers mapping
            single_reflection_tables.extend(table.split_by_experiment_id())
        else:
            single_reflection_tables.append(table)
    return single_reflection_tables


def renumber_table_id_columns(reflection_tables):
    """Renumber the id columns in the tables from 0..n-1

    If set, the experiment identifiers mapping is updated.

    Args:
        reflection_tables (list): A list of reflection tables

    Returns:
        (list): A list of reflection tables.
    """
    new_id_ = 0
    for table in reflection_tables:
        table_id_values = sorted(list(set(table["id"]).difference({-1})))
        highest_new_id = new_id_ + len(table_id_values) - 1
        expt_ids_dict = table.experiment_identifiers()
        new_ids_dict = {}
        new_id_ = highest_new_id
        while table_id_values:
            val = table_id_values.pop()
            sel = table["id"] == val
            if val in expt_ids_dict:
                # only delete here, add new at end to avoid clashes of new/old ids
                new_ids_dict[new_id_] = expt_ids_dict[val]
                del expt_ids_dict[val]
            table["id"].set_selected(sel.iselection(), new_id_)
            new_id_ -= 1
        new_id_ = highest_new_id + 1
        if new_ids_dict:
            for i, v in new_ids_dict.items():
                expt_ids_dict[i] = v
    return reflection_tables


def parse_multiple_datasets(reflection_tables):
    single_reflection_tables = split_reflection_tables_on_ids(reflection_tables)
    return renumber_table_id_columns(single_reflection_tables)


def sort_tables_to_experiments_order(reflection_tables, experiments):
    """If experiment identifiers are set, sort the order of reflection tables
    to match the order of the experiments.

    example for several single datasets in order
    input [r1("0"), r2("1")], [exp("0"), exp("1")]
    returns [r1("0"), r2("1")]

    example for several single datasets out of order
    input [r1("1"), r2("0")], [exp("0"), exp("1")]
    returns [r2("0"), r1("1")]

    example including a multi-dataset reflection table
    (e.g. datasets from command line: d1&2.refl d0.refl d0.expt d1&2.expt)
    input [r1("1", "2"), r2("0")], [exp("0"), exp("1"), exp("2")]
    returns [r2("0"), r1("1", "2")]

    Args:
        reflection_tables (list): A list of reflection tables (may contain multiple
            datasets per table)
        experiments: An ExperimentList

    Returns:
        (list): A list of sorted reflection tables to match the experiments order.
    """
    identifiers_list = []
    identifiers_by_table_idx = {}
    exp_id_to_table_idx = {}
    for i, table in enumerate(reflection_tables):
        id_values = table.experiment_identifiers().values()
        if id_values:
            identifiers_list.extend(id_values)
            identifiers_by_table_idx[i] = id_values
            for id_ in id_values:
                exp_id_to_table_idx[id_] = i

    expt_identiers = list(experiments.identifiers())
    if len(identifiers_list) == len(experiments) and set(identifiers_list) == set(
        expt_identiers
    ):  # all set so can now do the rearrangement
        sorted_tables = []
        for id_ in expt_identiers:
            # check if still need to add a table for this id
            if id_ in exp_id_to_table_idx:
                table_idx = exp_id_to_table_idx[id_]  # find table index
                sorted_tables.append(reflection_tables[table_idx])
                # now delete remaining ids from id>table map so that don't add twice
                for exp_id in identifiers_by_table_idx[table_idx]:
                    del exp_id_to_table_idx[exp_id]
        assert len(sorted_tables) == len(reflection_tables)
        return sorted_tables
    # else, identifiers not fully set, just return the tables.
    return reflection_tables


def assign_unique_identifiers(experiments, reflections, identifiers=None):
    """
    Assign unique experiment identifiers to experiments and reflections lists.

    If experiment identifiers are not set for some datasets, then new unique
    identifiers are given to those. The id column values are preserved.

    Args:
        experiments: An ExperimentList
        reflections (list): A list of reflection tables

    Returns:
        (tuple): tuple containing:
            experiments: The updated ExperimentList
            reflections (list): A list of the updated reflection tables

    Raises:
        ValueError: If the number of reflection tables and experiments are
            unequal (and identifiers if specified). Also raised if the existing
            identifiers are corrupted.
    """
    if len(experiments) != len(reflections):
        raise ValueError(
            "The experiments and reflections lists are unequal in length: %s & %s"
            % (len(experiments), len(reflections))
        )
    for i, table in enumerate(reflections):
        n_datasets = len(set(table["id"]).difference({-1}))
        if n_datasets > 1:
            raise ValueError(
                "Reflection table %s contains %s datasets (must contain a single dataset)"
                % (i, n_datasest)
            )
    # if identifiers given, use these to set the identifiers
    if identifiers:
        if len(identifiers) != len(reflections):
            raise ValueError(
                "The identifiers and reflections lists are unequal in length: %s & %s"
                % (len(identifiers), len(reflections))
            )
        for exp, refl, identifier in zip(experiments, reflections, identifiers):
            exp.identifier = identifier
            id_ = list(set(refl["id"]).difference({-1}))[
                0
            ]  # earlier check ensures only one
            refl.experiment_identifiers()[id_] = identifier
            refl.clean_experiment_identifiers_map()
    # Validate the existing identifiers, or the ones just set
    used_str_ids = []
    for exp, refl in zip(experiments, reflections):
        if exp.identifier != "":
            if list(refl.experiment_identifiers().values()) != [exp.identifier]:
                raise ValueError(
                    "Corrupted identifiers, please check input: in reflections: %s, in experiment: %s"
                    % (list(refl.experiment_identifiers().values()), exp.identifier)
                )
            used_str_ids.append(exp.identifier)

    if len(set(used_str_ids)) != len(reflections):
        # Generate a uuid for any identifiers not set.
        for exp, refl in zip(experiments, reflections):
            if exp.identifier == "":
                strid = str(uuid.uuid4())
                exp.identifier = strid
                id_ = list(set(refl["id"]).difference({-1}))[0]
                refl.experiment_identifiers()[id_] = strid
                refl.clean_experiment_identifiers_map()
    return experiments, reflections


def select_datasets_on_ids(
    experiments, reflection_table_list, exclude_datasets=None, use_datasets=None
):
    """
    Select a subset of experiments and reflection tables based on identifiers.

    This performs a similar function to the select/remove_on_experiment_identifiers
    methods of ExperimentList and reflection_table, with additional logic to handle
    the case of a list of reflection tables, rather than a single one and to catch
    bad input. Does not require reflection tables containing data from multiple
    experiments to be split.

    Args:
        experiments: An ExperimentList
        reflection_table_list (list): a list of reflection tables
        exclude_datasets (list): a list of experiment_identifiers to exclude
        use_datasets (list): a list of experiment_identifiers to use

    Returns:
        (tuple): tuple containing:
            experiments: The updated ExperimentList
            list_of_reflections (list): A list of the updated reflection tables

    Raises:
        ValueError: If both use_datasets and exclude datasets are used, if not all
            experiment identifiers are set, if an identifier in exclude_datasets
            or use_datasets is not in the list.
    """
    if not use_datasets and not exclude_datasets:
        return experiments, reflection_table_list
    if use_datasets and exclude_datasets:
        raise ValueError(
            "The options use_datasets and exclude_datasets cannot be used in conjuction."
        )
    if experiments.identifiers().count("") > 0:
        raise ValueError(
            """
            Not all experiment identifiers set in the ExperimentList.
            Current identifiers set as: %s"""
            % list(experiments.identifiers())
        )
    list_of_reflections = []
    if use_datasets:
        total_found = 0
        for reflection_table in reflection_table_list:
            expids = list(reflection_table.experiment_identifiers().values())
            expids_in_this_table = list(set(use_datasets).intersection(set(expids)))
            total_found += len(expids_in_this_table)
            if expids_in_this_table:  # if none, then no datasets wanted from table
                list_of_reflections.append(
                    reflection_table.select_on_experiment_identifiers(
                        expids_in_this_table
                    )
                )
        if total_found != len(use_datasets):
            raise ValueError(
                """Attempting to select datasets based on identifiers that
are not found in the experiment list / reflection tables."""
            )
        experiments.select_on_experiment_identifiers(use_datasets)
    elif exclude_datasets:
        if len(
            set(exclude_datasets).intersection(set(experiments.identifiers()))
        ) != len(set(exclude_datasets)):
            raise ValueError(
                """Attempting to exclude datasets based on identifiers that
are not found in the experiment list / reflection tables."""
            )
        for reflection_table in reflection_table_list:
            expids = list(reflection_table.experiment_identifiers().values())
            expids_in_this_table = list(set(exclude_datasets).intersection(set(expids)))
            if expids_in_this_table:  # only append if data left after removing
                r_table = reflection_table.remove_on_experiment_identifiers(
                    expids_in_this_table
                )
                if r_table.size() > 0:
                    list_of_reflections.append(r_table)
            else:
                list_of_reflections.append(reflection_table)
        experiments.remove_on_experiment_identifiers(exclude_datasets)
    return experiments, list_of_reflections

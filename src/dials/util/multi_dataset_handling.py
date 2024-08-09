"""
Module of functions for handling operations on lists of reflection tables
and experiment lists.
"""


from __future__ import annotations

import copy
import logging

from orderedset import OrderedSet

from dxtbx.model import ExperimentList
from dxtbx.util import ersatz_uuid4

from dials.array_family import flex

logger = logging.getLogger("dials")

import iotbx.phil

phil_scope = iotbx.phil.parse(
    """
  dataset_selection {
    use_datasets = None
      .type = ints
      .help = "Choose a subset of datasets, based on the dataset id (as defined
               in the reflection table), to use from a multi-dataset input."
      .expert_level = 2
    include_datasets = None
      .type = ints
      .help = "An alias for use_datasets option"
      .expert_level = 3
    exclude_datasets = None
      .type = ints
      .help = "Choose a subset of datasets, based on the dataset id (as defined
               in the reflection table), to exclude from a multi-dataset input."
      .expert_level = 2
  }
"""
)


class Expeditor(object):
    def __init__(self, experiments, reflection_tables=None):

        self.experiments = experiments
        self.reflection_table = None

        self.reflections_with_id_minus_1 = flex.reflection_table()
        if reflection_tables:
            if len(reflection_tables) > 1:
                n_imagesets = 0
                for t in reflection_tables:
                    if "imageset_id" in t:
                        t["imageset_id"] += n_imagesets
                        n_imagesets += len(set(t["imageset_id"]))
                self.reflection_table = flex.reflection_table.concat(reflection_tables)
            else:
                self.reflection_table = reflection_tables[0]
            if -1 in set(self.reflection_table["id"]):
                sel = self.reflection_table["id"] == -1
                self.reflections_with_id_minus_1 = self.reflection_table.select(sel)
                self.reflection_table.del_selected(sel)
        self.crystal_locs = [
            i for i, expt in enumerate(self.experiments) if expt.crystal
        ]
        self.input_has_crystalless_expts = bool(
            len(self.crystal_locs) != len(self.experiments)
        )
        self.crystalless_reflection_tables = []

    def get_unindexed(self):
        return self.reflections_with_id_minus_1

    def filter_experiments_with_crystals(self):
        if not self.input_has_crystalless_expts:
            if not self.reflection_table:
                return self.experiments, None
            tables = parse_multiple_datasets([self.reflection_table])
            if not self.experiments.identifiers():
                assign_unique_identifiers(tables, self.experiments)
            return self.experiments, tables
        if not self.experiments.identifiers():
            for i, expt in enumerate(self.experiments):
                strid = ersatz_uuid4()
                expt.identifier = strid
                if self.reflection_table:
                    self.reflection_table.experiment_identifiers()[i] = strid

        expts_with_crystals = ExperimentList(
            [expt for expt in self.experiments if expt.crystal]
        )
        reflection_table = None
        if self.reflection_table:
            reflection_table = self.reflection_table.select_on_experiment_identifiers(
                expts_with_crystals.identifiers()
            )
            reflection_table.reset_ids()

            # record some quantities to help later when recombining
            self.crystalless_reflection_tables = [
                self.reflection_table.select_on_experiment_identifiers(
                    [expt.identifier]
                )
                if not expt.crystal
                else None
                for expt in self.experiments
            ]
            reflection_table = reflection_table.split_by_experiment_id()
        return expts_with_crystals, reflection_table

    def combine_experiments_for_output(
        self, experiments, reflection_tables=None, include_unindexed=False
    ):
        # When reflection tables exist, handle two types of input.
        # Either the program returns a list of reflections tables, one for each experiment,
        # or a combined table for all experiments.
        if reflection_tables:
            if len(experiments) != len(reflection_tables):
                if len(reflection_tables) == 1:
                    reflection_tables[0].assert_experiment_identifiers_are_consistent(
                        experiments
                    )
                    unassigned = None
                    new_unassigned_sel = reflection_tables[0]["id"] == -1
                    if any(new_unassigned_sel):
                        logger.info(f"{new_unassigned_sel.count(True)} unassigned")
                        unassigned = reflection_tables[0].select(new_unassigned_sel)
                    reflection_tables = reflection_tables[0].split_by_experiment_id()
                    if unassigned and include_unindexed:
                        reflection_tables.append(unassigned)
                else:
                    raise ValueError("Mismatch")
        imagesets = self.experiments.imagesets()
        if not self.input_has_crystalless_expts:
            if reflection_tables:
                for expt, table in zip(experiments, reflection_tables):
                    if (
                        "imageset_id" in table and imagesets
                    ):  # strong check to guarantee
                        table["imageset_id"] = flex.int(
                            table.size(), imagesets.index(expt.imageset)
                        )
                if self.reflections_with_id_minus_1 and include_unindexed:
                    reflection_tables.append(self.reflections_with_id_minus_1)
                return experiments, flex.reflection_table.concat(reflection_tables)
            return experiments, None

        for i, expt, table in zip(self.crystal_locs, experiments, reflection_tables):
            other_expts_sharing_scan = self.experiments.where(scan=expt.scan)
            for j in other_expts_sharing_scan:
                self.experiments[j].beam = expt.beam
                self.experiments[j].detector = expt.detector
                self.experiments[j].goniometer = expt.goniometer
            self.experiments[i] = expt
            self.crystalless_reflection_tables[i] = table

        for expt, table in zip(self.experiments, self.crystalless_reflection_tables):
            if "imageset_id" in table:
                table["imageset_id"] = flex.int(
                    table.size(), imagesets.index(expt.imageset)
                )

        reflections = flex.reflection_table.concat(self.crystalless_reflection_tables)
        reflections.assert_experiment_identifiers_are_consistent(self.experiments)
        return self.experiments, reflections

    def generate_experiments_with_updated_crystal(self, experiments, crystal_id):
        # special function to handle the awkward case of dials.rbs
        crystals = [
            expt.crystal for expt in self.experiments if expt.crystal is not None
        ]
        crystal = crystals[crystal_id]
        expts_ids = self.experiments.where(crystal=crystal)
        # make a copy
        elist = copy.copy(self.experiments)

        for id_, new_expt in zip(expts_ids, experiments):
            elist[id_].crystal = new_expt.crystal
            # also copy across beam, detector and gonio for anything sharing those models
            scan = elist[id_].scan
            j_expt = elist.where(scan=scan)
            for j in j_expt:
                elist[j].beam = new_expt.beam
                elist[j].goniometer = new_expt.goniometer
                elist[j].detector = new_expt.detector
        return elist


def generate_experiment_identifiers(experiments, identifier_type="uuid"):
    """Generate unique identifiers for each experiment."""
    if identifier_type == "uuid":
        for expt in experiments:
            expt.identifier = ersatz_uuid4()
    elif identifier_type == "timestamp":
        pass


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
        all_blanks = all(len(id_value) == 0 for id_value in id_values)
        if id_values and not all_blanks:
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
        if not table:
            continue
        table_id_values = sorted(set(table["id"]).difference({-1}), reverse=True)
        highest_new_id = new_id_ + len(table_id_values) - 1
        expt_ids_dict = table.experiment_identifiers()
        new_ids_dict = {}
        new_id_ = highest_new_id
        orig_id = copy.deepcopy(table["id"])
        for val in table_id_values:
            sel = orig_id == val
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


def parse_multiple_datasets(reflections):
    """
    Split a list of multi-dataset reflection tables, selecting on id

    If duplicate id values are found, the id columns are renumbered from 0..n-1,
    taking care of experiment identifiers if these are set.

    Args:
        reflections (list): a list of reflection tables, each of which may contain
            multiple datasets

    Returns:
        (list): a list of reflection tables corresponding to single datasets
    """
    single_reflection_tables = []
    dataset_id_list = []
    for refl_table in reflections:
        dataset_ids = set(refl_table["id"]).difference({-1})
        dataset_id_list.extend(list(dataset_ids))
        if len(dataset_ids) > 1:
            logger.info(
                "Detected existence of a multi-dataset reflection table \n"
                "containing %s datasets. \n",
                len(dataset_ids),
            )
            # FIXME fix split_by_experiment_id so that don't need to filter
            # unindxeded reflections here to get rid of id = -1
            if -1 in refl_table["id"]:
                refl_table = refl_table.select(refl_table["id"] != -1)
            result = refl_table.split_by_experiment_id()
            single_reflection_tables.extend(result)
        else:
            single_reflection_tables.append(refl_table)
    if len(dataset_id_list) != len(set(dataset_id_list)):  # need to reset some ids
        logger.warning(
            "Duplicate dataset ids found in different reflection tables. "
            "These will be treated as coming from separate datasets, and "
            "new dataset ids will be assigned for the whole dataset."
        )
        for new_id, (r, old_id) in enumerate(
            zip(single_reflection_tables, dataset_id_list)
        ):
            r["id"] = flex.int(r.size(), new_id)
            if list(r.experiment_identifiers()):  # if identifiers, need to update
                expid = r.experiment_identifiers()[old_id]
                del r.experiment_identifiers()[old_id]
                r.experiment_identifiers()[new_id] = expid
    return single_reflection_tables


def assign_unique_identifiers(experiments, reflections, identifiers=None):
    """
    Assign unique experiment identifiers to experiments and reflections lists.

    If experiment identifiers are not set for some datasets, then new unique
    identifiers are given to those, and the 'id' column for all reflection tables
    are set sequentially from 0..n-1.

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
    # if identifiers given, use these to set the identifiers
    if identifiers:
        if len(identifiers) != len(reflections):
            raise ValueError(
                "The identifiers and reflections lists are unequal in length: %s & %s"
                % (len(identifiers), len(reflections))
            )
        for i, (exp, refl) in enumerate(zip(experiments, reflections)):
            exp.identifier = identifiers[i]
            for k in refl.experiment_identifiers().keys():
                del refl.experiment_identifiers()[k]
            refl.experiment_identifiers()[i] = identifiers[i]
            refl["id"] = flex.int(refl.size(), i)
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
        # if not all set, then need to fill in the rest. Keep the identifier if
        # it is already set, and reset table id column from 0..n-1
        for i, (exp, refl) in enumerate(zip(experiments, reflections)):
            if exp.identifier == "":
                strid = ersatz_uuid4()
                exp.identifier = strid
                refl.experiment_identifiers()[i] = strid
            else:
                k = list(refl.experiment_identifiers().keys())[0]
                expid = list(refl.experiment_identifiers().values())[0]
                del refl.experiment_identifiers()[k]
                refl.experiment_identifiers()[i] = expid
            refl["id"] = flex.int(refl.size(), i)
    return experiments, reflections


def select_datasets_on_ids(
    experiments, reflection_table_list, exclude_datasets=None, use_datasets=None
):
    # transform dataset ids to identifiers
    if not use_datasets and not exclude_datasets:
        return experiments, reflection_table_list
    if use_datasets and exclude_datasets:
        raise ValueError(
            "The options use_datasets and exclude_datasets cannot be used in conjunction."
        )

    id_map = {}
    for table in reflection_table_list:
        id_map.update(table.experiment_identifiers())
    if exclude_datasets:
        identifiers_to_exclude = []
        for k in exclude_datasets:
            if int(k) not in id_map:
                raise ValueError(
                    """Attempting to select datasets based on identifiers that
are not found in the experiment list / reflection tables."""
                )
            identifiers_to_exclude.append(id_map[int(k)])
        return select_datasets_on_identifiers(
            experiments, reflection_table_list, exclude_datasets=identifiers_to_exclude
        )
    elif use_datasets:
        identifiers_to_use = []
        for k in use_datasets:
            if int(k) not in id_map:
                raise ValueError(
                    """Attempting to select datasets based on identifiers that
are not found in the experiment list / reflection tables."""
                )
            identifiers_to_use.append(id_map[int(k)])
        return select_datasets_on_identifiers(
            experiments, reflection_table_list, use_datasets=identifiers_to_use
        )


def select_datasets_on_identifiers(
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
            "The options use_datasets and exclude_datasets cannot be used in conjunction."
        )
    if experiments.identifiers().count("") > 0:
        raise ValueError(
            f"""
            Not all experiment identifiers set in the ExperimentList.
            Current identifiers set as: {list(experiments.identifiers())}"""
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


def update_imageset_ids(experiments, reflections):
    """For a list of input experiments and reflections (each containing one
    sweep), update or add the imageset_id column to the data to match the order
    in the experiment list.

    This means that when the reflection tables are combined, the data is correct.
    """
    # input a list of ordered matching experiments and reflection tables.

    next_iset_id = 0
    imagesets_found = OrderedSet()
    for expt, table in zip(experiments, reflections):
        if "imageset_id" in table:
            assert len(set(table["imageset_id"])) == 1
        iset = expt.imageset
        if iset not in imagesets_found:
            imagesets_found.add(iset)
            table["imageset_id"] = flex.int(table.size(), next_iset_id)
            next_iset_id += 1
        else:
            iset_id = imagesets_found.index(iset)
            table["imageset_id"] = flex.int(table.size(), iset_id)
    return reflections

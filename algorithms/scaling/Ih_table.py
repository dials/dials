"""
A datastructure for summing over groups of symmetry equivalent reflections.

This module defines a blocked datastructures for summing over groups of
symmetry equivalent reflections, as required for scaling.
"""
from __future__ import absolute_import, division, print_function

from orderedset import OrderedSet
from dials.array_family import flex
from cctbx import miller, crystal, uctbx
from scitbx import sparse


def map_indices_to_asu(miller_indices, space_group):
    """Map the indices to the asymmetric unit."""
    crystal_symmetry = crystal.symmetry(space_group=space_group)
    miller_set = miller.set(
        crystal_symmetry=crystal_symmetry, indices=miller_indices, anomalous_flag=False
    )
    miller_set_in_asu = miller_set.map_to_asu()
    return miller_set_in_asu.indices()


def get_sorted_asu_indices(asu_indices, space_group):
    """Return the sorted asu indices and the permutation selection."""
    crystal_symmetry = crystal.symmetry(space_group=space_group)
    miller_set_in_asu = miller.set(
        crystal_symmetry=crystal_symmetry, indices=asu_indices, anomalous_flag=False
    )
    permuted = miller_set_in_asu.sort_permutation(by_value="packed_indices")
    sorted_asu_miller_index = asu_indices.select(permuted)
    return sorted_asu_miller_index, permuted


class IhTable(object):
    """
    A class to manage access to Ih_table blocks.

    The idea here is to split the data into blocks to allow parallelized
    computations, but within the blocks the data are sorted by dataset.
    In each block, there exists a block_selection_list which contains the indices
    for each dataset from the input reflection table.

    This class acts as a 'master' to setup the block structure and control access
    to the underlying blocks - only metadata is kept in this class after
    initialisation, the reflections etc are all contained in the blocks.
    To set data in the blocks, methods are provided by the master, e.g
    set_intensities(intensities, block_id) which are delegated down to the
    appropriate block.

    Attributes:
        space_group: The space group for the dataset.
        Ih_table_blocks (list): A list of IhTableBlock instances. All symmetry
            equivalent reflections are recorded in the same block, to allow
            splitting of the dataset for parallelized computations.
        nblocks (int): The number of blocks in the Ih_table_blocks list.
        blocked_selection_list (list): A list of lists. bsl[i][j] is the selection
            list for block i, dataset j.
        n_datasets: The number of input reflection tables used to make the Ih_table.
        size: The number of reflections across all blocks
        asu_index_dict (dict): A dictionary, key: asu_miller_index, value tuple
            containing group_id and block_id (where group id is the group index
            within its block).

    """

    id_ = "IhTable"

    def __init__(
        self,
        reflection_tables,
        space_group,
        indices_lists=None,
        nblocks=1,
        free_set_percentage=0,
        free_set_offset=0,
        additional_cols=None,
    ):
        """
        Distribute the input data into the required structure.

        The reflection data can be split into blocks, while the relevant
        metadata is also generated.

        A list of flex.size_t indices can be provided - this allows the
        reflection table data to maintain a reference to a dataset from which
        it was selecte; these will be used when making the block selections.
        e.g selection = flex.bool([True, False, True])
            r_1 = r_master.select(selection)
            indices_list = selection.iselection() = flex.size_t([0, 2])
            then the block selection will contain 0 and 2 to refer back
            to the location of the data in r_master.
        """
        if indices_lists:
            assert len(indices_lists) == len(reflection_tables)
        self.asu_index_dict = {}
        self.space_group = space_group
        self.n_work_blocks = nblocks
        self.n_datasets = len(reflection_tables)
        self.Ih_table_blocks = []
        self.blocked_selection_list = []
        self.properties_dict = {
            "n_unique_in_each_block": [],
            "n_reflections_in_each_block": {},
            "miller_index_boundaries": [],
        }
        self._determine_required_block_structures(reflection_tables, self.n_work_blocks)
        self._create_empty_Ih_table_blocks()
        for i, table in enumerate(reflection_tables):
            if indices_lists:
                self._add_dataset_to_blocks(
                    i, table, indices_lists[i], additional_cols=additional_cols
                )
            else:
                self._add_dataset_to_blocks(i, table, additional_cols=additional_cols)
        self.generate_block_selections()
        self.free_Ih_table = None
        if free_set_percentage > 0:
            self.extract_free_set(free_set_percentage, offset=free_set_offset)
            self.free_Ih_table = True
        self.calc_Ih()

    def update_data_in_blocks(self, data, dataset_id, column="intensity"):
        """
        Update a given column across all blocks for a given dataset.

        Given an array of data (of the same size as the input reflection
        table) and the name of the column, use the internal data to split
        this up and set in individual blocks.
        """
        assert column in ["intensity", "variance", "inverse_scale_factor"]
        assert dataset_id in range(0, self.n_datasets)
        # split up data for blocks
        for block in self.blocked_data_list:
            data_for_block = data.select(block.block_selections[dataset_id])
            start = block.dataset_info[dataset_id]["start_index"]
            end = block.dataset_info[dataset_id]["end_index"]
            sel = flex.size_t(range(start, end))
            block.Ih_table[column].set_selected(sel, data_for_block)
            if column == "variance":
                block.Ih_table["weights"].set_selected(sel, 1.0 / data_for_block)

    def get_block_selections_for_dataset(self, dataset):
        """Generate the block selection list for a given dataset."""
        assert dataset in range(self.n_datasets)
        if self.free_Ih_table:
            return [
                self.blocked_selection_list[i][dataset]
                for i in range(self.n_work_blocks + 1)
            ]
        return [
            self.blocked_selection_list[i][dataset] for i in range(self.n_work_blocks)
        ]

    @property
    def size(self):
        """Sum the sizes of all work blocks to give the total number of reflections."""
        if self.free_Ih_table:
            return sum([block.size for block in self.Ih_table_blocks[:-1]])
        return sum([block.size for block in self.Ih_table_blocks])

    def generate_block_selections(self):
        """Generate and set an updated blocked_selection_list."""
        self.blocked_selection_list = [
            block.block_selections for block in self.Ih_table_blocks
        ]

    def update_weights(self, block_id=None):
        """Update weights (to allow iterative updating, not implemented)."""

    def update_error_model(self, error_model):
        """Update the error model in the blocks."""
        for block in self.Ih_table_blocks:
            block.update_error_model(error_model)

    def reset_error_model(self):
        """Reset the weights in the blocks."""
        for block in self.Ih_table_blocks:
            block.reset_error_model()

    @property
    def blocked_data_list(self):
        """Return the list of IhTableBlock instances."""
        return self.Ih_table_blocks

    def set_intensities(self, intensities, block_id):
        """Set the intensities for a given block."""
        self.Ih_table_blocks[block_id].Ih_table["intensity"] = intensities

    def set_derivatives(self, derivatives, block_id):
        """Set the derivatives matrix for a given block."""
        self.Ih_table_blocks[block_id].derivatives = derivatives

    def set_inverse_scale_factors(self, new_scales, block_id):
        """Set the inverse scale factors for a given block."""
        self.Ih_table_blocks[block_id].inverse_scale_factors = new_scales

    def set_variances(self, new_variances, block_id):
        """Set the variances and weights for a given block."""
        self.Ih_table_blocks[block_id].Ih_table["variance"] = new_variances
        self.Ih_table_blocks[block_id].Ih_table["weights"] = 1.0 / new_variances

    def calc_Ih(self, block_id=None):
        """Calculate the latest value of Ih, for a given block or for all blocks."""
        if block_id is not None:
            self.Ih_table_blocks[block_id].calc_Ih()
        else:
            for block in self.Ih_table_blocks:
                block.calc_Ih()

    def _determine_required_block_structures(self, reflection_tables, nblocks=1):
        """
        Inspect the input to determine how to split into blocks.

        Extract the asu miller indices from the reflection table and
        add data to the asu_index_dict and properties dict.
        """
        joint_asu_indices = flex.miller_index()
        for table in reflection_tables:
            if not "asu_miller_index" in table:
                table["asu_miller_index"] = map_indices_to_asu(
                    table["miller_index"], self.space_group
                )
            joint_asu_indices.extend(table["asu_miller_index"])
        sorted_joint_asu_indices, _ = get_sorted_asu_indices(
            joint_asu_indices, self.space_group
        )

        asu_index_set = OrderedSet(sorted_joint_asu_indices)
        n_unique_groups = len(asu_index_set)
        # also record how many unique groups go into each block
        group_boundaries = [int(i * n_unique_groups / nblocks) for i in range(nblocks)]
        group_boundaries.append(n_unique_groups)

        next_boundary = group_boundaries[1]
        block_id = 0
        group_id_in_block_i = 0
        for i, index in enumerate(asu_index_set):
            if i == next_boundary:
                self.properties_dict["n_unique_in_each_block"].append(
                    group_id_in_block_i
                )
                self.properties_dict["miller_index_boundaries"].append(index)
                block_id += 1
                next_boundary = group_boundaries[block_id + 1]
                group_id_in_block_i = 0
            self.asu_index_dict[index] = (group_id_in_block_i, block_id)
            group_id_in_block_i += 1
        # record the number in the last block
        self.properties_dict["n_unique_in_each_block"].append(group_id_in_block_i)
        self.properties_dict["miller_index_boundaries"].append((10000, 10000, 10000))
        # ^ to avoid bounds checking when in last group
        # need to know how many reflections will be in each block also

        block_id = 0
        idx_prev = 0
        boundary = self.properties_dict["miller_index_boundaries"][0]
        for i, index in enumerate(sorted_joint_asu_indices):
            if index == boundary:
                n_in_prev_group = i - idx_prev
                self.properties_dict["n_reflections_in_each_block"][
                    block_id
                ] = n_in_prev_group
                block_id += 1
                boundary = self.properties_dict["miller_index_boundaries"][block_id]
                idx_prev = i
        self.properties_dict["n_reflections_in_each_block"][block_id] = (
            len(sorted_joint_asu_indices) - idx_prev
        )

    def _create_empty_Ih_table_blocks(self):
        for n in range(self.n_work_blocks):
            n_refl_in_block = self.properties_dict["n_reflections_in_each_block"][n]
            n_groups_in_block = self.properties_dict["n_unique_in_each_block"][n]
            self.Ih_table_blocks.append(
                IhTableBlock(
                    n_groups=n_groups_in_block,
                    n_refl=n_refl_in_block,
                    n_datasets=self.n_datasets,
                )
            )

    def _add_dataset_to_blocks(
        self, dataset_id, reflections, indices_array=None, additional_cols=None
    ):
        sorted_asu_indices, perm = get_sorted_asu_indices(
            reflections["asu_miller_index"], self.space_group
        )
        r = flex.reflection_table()
        r["intensity"] = reflections["intensity"]
        r["asu_miller_index"] = reflections["asu_miller_index"]
        r["variance"] = reflections["variance"]
        r["inverse_scale_factor"] = reflections["inverse_scale_factor"]
        if isinstance(additional_cols, list):
            for col in additional_cols:
                if col in reflections:
                    r[col] = reflections[col]
        if indices_array:
            r["loc_indices"] = indices_array
        else:
            r["loc_indices"] = flex.size_t(range(r.size()))
        r = r.select(perm)
        r["dataset_id"] = flex.int(r.size(), dataset_id)
        # if data are sorted by asu_index, then up until boundary, should be in same
        # block (still need to read group_id though)

        # sort data, get group ids and block_ids
        group_ids = flex.int([])
        boundary = self.properties_dict["miller_index_boundaries"][0]
        boundary_id = 0
        boundaries_for_this_datset = [0]  # use to slice
        # make this a c++ method for speed?
        for i, index in enumerate(sorted_asu_indices):
            while index >= boundary:
                boundaries_for_this_datset.append(i)
                boundary_id += 1
                boundary = self.properties_dict["miller_index_boundaries"][boundary_id]
            group_id, _ = self.asu_index_dict[index]
            group_ids.append(group_id)
        while len(boundaries_for_this_datset) < self.n_work_blocks + 1:
            # catch case where last boundaries aren't reached
            boundaries_for_this_datset.append(len(sorted_asu_indices))
        # so now have group ids as well for individual dataset
        if self.n_work_blocks == 1:
            self.Ih_table_blocks[0].add_data(dataset_id, group_ids, r)
        else:
            for i, val in enumerate(boundaries_for_this_datset[:-1]):
                start = val
                end = boundaries_for_this_datset[i + 1]
                self.Ih_table_blocks[i].add_data(
                    dataset_id, group_ids[start:end], r[start:end]
                )

    def extract_free_set(self, free_set_percentage, offset=0):
        """Extract a free set from all blocks."""
        assert not self.free_Ih_table
        interval_between_groups = int(100 / free_set_percentage)
        free_reflection_table = flex.reflection_table()
        free_indices = flex.size_t()
        for j, block in enumerate(self.Ih_table_blocks):
            n_groups = block.h_index_matrix.n_cols
            groups_for_free_set = flex.bool(n_groups, False)
            for_free = flex.size_t(
                [i for i in range(0 + offset, n_groups, interval_between_groups)]
            )
            groups_for_free_set.set_selected(for_free, True)
            free_block = block.select_on_groups(groups_for_free_set)
            free_reflection_table.extend(free_block.Ih_table)
            for sel in free_block.block_selections:
                free_indices.extend(sel)
            self.Ih_table_blocks[j] = block.select_on_groups(~groups_for_free_set)
            # Now need to update dataset_info dict.
            removed_from_each_dataset = [
                (free_block.Ih_table["dataset_id"] == i).count(True)
                for i in range(0, block.n_datasets)
            ]
            n_removed = 0
            for i in range(0, self.Ih_table_blocks[j].n_datasets):
                self.Ih_table_blocks[j].dataset_info[i]["start_index"] -= n_removed
                n_removed += removed_from_each_dataset[i]
                self.Ih_table_blocks[j].dataset_info[i]["end_index"] -= n_removed
        self.blocked_selection_list = [
            block.block_selections for block in self.Ih_table_blocks
        ]
        # now split by dataset and use to instantiate another Ih_table
        datasets = set(free_reflection_table["dataset_id"])
        tables = []
        indices_lists = []
        for id_ in datasets:
            dataset_sel = free_reflection_table["dataset_id"] == id_
            tables.append(free_reflection_table.select(dataset_sel))
            indices_lists.append(free_indices.select(dataset_sel))
        free_Ih_table = IhTable(tables, self.space_group, indices_lists, nblocks=1)
        # add to blocks list and selection list
        self.Ih_table_blocks.append(free_Ih_table.blocked_data_list[0])
        self.blocked_selection_list.append(free_Ih_table.blocked_selection_list[0])

    def as_miller_array(self, unit_cell, return_free_set_data=False):
        """Get a scaled miller array from the Ih_table and an experiment."""
        blocked_data_list = self.blocked_data_list

        if self.free_Ih_table:
            if return_free_set_data:
                blocked_data_list = [blocked_data_list[-1]]
            else:
                blocked_data_list = blocked_data_list[:-1]
        if len(blocked_data_list) > 1:
            joint_table = flex.reflection_table()
            for block in blocked_data_list:
                # better to just create many miller arrays and join them?
                refl_for_joint_table = flex.reflection_table()
                for col in [
                    "asu_miller_index",
                    "intensity",
                    "inverse_scale_factor",
                    "variance",
                ]:
                    refl_for_joint_table[col] = block.Ih_table[col]
                joint_table.extend(refl_for_joint_table)
        else:
            joint_table = blocked_data_list[0].Ih_table
        # Filter out negative scale factors to avoid merging statistics errors.
        return _reflection_table_to_iobs(joint_table, unit_cell, self.space_group)


class IhTableBlock(object):
    """
    A datastructure for efficient summations over symmetry equivalent reflections.

    This contains a reflection table, sorted by dataset, called the Ih_table,
    a h_index_matrix (sparse) for efficiently calculating sums over symmetry
    equivalent reflections as well as 'block_selections' which relate the order
    of the data to the initial reflection tables used to initialise the (master)
    IhTable.

    Attributes:
        Ih_table: A reflection table, containing I, g, w, var, Ih,
            asu_miller_index, loc_indices and dataset_id.
        block_selections: A list of flex.size_t arrays of indices, that can be
            used to select and reorder data from the input reflection tables to
            match the order in the Ih_table.
        h_index_matrix: A sparse matrix used to sum over groups of equivalent
            reflections by multiplication. Sum_h I = I * h_index_matrix. The
            dimension is n_refl by n_groups; each row has a single nonzero
            entry with a value of 1.
        h_expand_matrix: The transpose of the h_index_matrix, used to expand an
            array of values for symmetry groups into an array of size n_refl.
        derivatives: A matrix of derivatives of the reflections wrt the model
            parameters.

    """

    def __init__(self, n_groups, n_refl, n_datasets=1):
        """Create empty datastructures to which data can later be added."""
        self.Ih_table = flex.reflection_table()
        self.block_selections = [None] * n_datasets
        self.h_index_matrix = sparse.matrix(n_refl, n_groups)
        self._setup_info = {"next_row": 0, "next_dataset": 0, "setup_complete": False}
        self.dataset_info = {}
        self.n_datasets = n_datasets
        self.h_expand_matrix = None
        self.derivatives = None
        self.binner = None

    def add_data(self, dataset_id, group_ids, reflections):
        """
        Add data to all blocks for a given dataset.

        Add data to the Ih_table, write data to the h_index_matrix and
        add the loc indices to the block_selections list.
        """
        assert not self._setup_info[
            "setup_complete"
        ], """
No further data can be added to the IhTableBlock as setup marked complete."""
        assert (
            self._setup_info["next_row"] + len(group_ids) <= self.h_index_matrix.n_rows
        ), """
Not enough space left to add this data, please check for correct block initialisation."""
        assert (
            dataset_id == self._setup_info["next_dataset"]
        ), """
Datasets must be added in correct order: expected: %s, this dataset: %s""" % (
            self._setup_info["next_dataset"],
            dataset_id,
        )
        assert "asu_miller_index" in reflections
        for i, id_ in enumerate(group_ids):
            self.h_index_matrix[i + self._setup_info["next_row"], id_] = 1.0
        self.dataset_info[dataset_id] = {"start_index": self._setup_info["next_row"]}
        self._setup_info["next_row"] += len(group_ids)
        self._setup_info["next_dataset"] += 1
        self.dataset_info[dataset_id]["end_index"] = self._setup_info["next_row"]
        self.Ih_table.extend(reflections)
        if "loc_indices" in reflections:
            self.block_selections[dataset_id] = reflections["loc_indices"]
        else:
            self.block_selections[dataset_id] = flex.int(range(len(reflections)))
        if self._setup_info["next_dataset"] == len(self.block_selections):
            self._complete_setup()

    def _complete_setup(self):
        """Finish the setup of the Ih_table once all data has been added."""
        self.h_index_matrix.compact()
        assert (
            self._setup_info["next_row"] == self.h_index_matrix.n_rows
        ), """
Not all rows of h_index_matrix appear to be filled in IhTableBlock setup."""
        self.h_expand_matrix = self.h_index_matrix.transpose()
        self.Ih_table["weights"] = 1.0 / self.Ih_table["variance"]
        self._setup_info["setup_complete"] = True

    def select(self, sel):
        """Select a subset of the data, returning a new IhTableBlock object."""
        Ih_table = self.Ih_table.select(sel)
        h_idx_sel = self.h_expand_matrix.select_columns(sel.iselection())
        reduced_h_idx = h_idx_sel.transpose()
        unity = flex.double(reduced_h_idx.n_rows, 1.0)
        nz_col_sel = (unity * reduced_h_idx) > 0
        h_index_matrix = reduced_h_idx.select_columns(nz_col_sel.iselection())
        h_expand = h_index_matrix.transpose()
        newtable = IhTableBlock(n_groups=0, n_refl=0, n_datasets=self.n_datasets)
        newtable.Ih_table = Ih_table
        newtable.h_expand_matrix = h_expand
        newtable.h_index_matrix = h_index_matrix
        newtable.block_selections = []
        offset = 0
        for i in range(newtable.n_datasets):
            newtable.dataset_info[i] = {"start_index": offset}
            block_sel_i = self.block_selections[i]
            n_in_dataset_i = len(block_sel_i)
            newtable.block_selections.append(
                block_sel_i.select(sel[offset : offset + n_in_dataset_i])
            )
            offset += n_in_dataset_i
            newtable.dataset_info[i]["end_index"] = offset
        return newtable

    def select_on_groups(self, sel):
        """Select a subset of the unique groups, returning a new IhTableBlock."""
        reduced_h_idx = self.h_index_matrix.select_columns(sel.iselection())
        unity = flex.double(reduced_h_idx.n_cols, 1.0)
        nz_row_sel = (unity * reduced_h_idx.transpose()) > 0
        return self.select(nz_row_sel)

    def select_on_groups_isel(self, isel):
        """Select a subset of the unique groups, returning a new IhTableBlock."""
        reduced_h_idx = self.h_index_matrix.select_columns(isel)
        unity = flex.double(reduced_h_idx.n_cols, 1.0)
        nz_row_sel = (unity * reduced_h_idx.transpose()) > 0
        return self.select(nz_row_sel)

    def calc_Ih(self):
        """Calculate the current best estimate for Ih for each reflection group."""
        scale_factors = self.Ih_table["inverse_scale_factor"]
        gsq = (scale_factors ** 2) * self.Ih_table["weights"]
        sumgsq = gsq * self.h_index_matrix
        gI = (scale_factors * self.Ih_table["intensity"]) * self.Ih_table["weights"]
        sumgI = gI * self.h_index_matrix
        Ih = sumgI / sumgsq
        self.Ih_table["Ih_values"] = Ih * self.h_expand_matrix

    def update_error_model(self, error_model):
        """Update the scaling weights based on an error model."""
        sigmaprimesq = error_model.update_variances(
            self.Ih_table["variance"], self.Ih_table["intensity"]
        )
        self.Ih_table["weights"] = 1.0 / sigmaprimesq

    def reset_error_model(self):
        """Reset the weights to their initial value."""
        self.Ih_table["weights"] = 1.0 / self.Ih_table["variance"]

    def calc_nh(self):
        """Calculate the number of refls in the group to which the reflection belongs.

        This is a vector of length n_refl."""
        return (
            flex.double(self.size, 1.0) * self.h_index_matrix
        ) * self.h_expand_matrix

    def match_Ih_values_to_target(self, target_Ih_table):
        """
        Use an Ih_table as a target to set Ih values in this table.

        Given an Ih table as a target, the common reflections across the tables
        are determined and the Ih_values are set to those of the target. If no
        matching reflection is found, then the values are removed from the table.
        """
        assert target_Ih_table.n_work_blocks == 1
        target_asu_Ih_dict = dict(
            zip(
                target_Ih_table.blocked_data_list[0].asu_miller_index,
                target_Ih_table.blocked_data_list[0].Ih_values,
            )
        )
        new_Ih_values = flex.double(self.size, 0.0)
        location_in_unscaled_array = 0
        sorted_asu_indices, permuted = get_sorted_asu_indices(
            self.Ih_table["asu_miller_index"], target_Ih_table.space_group
        )
        for j, miller_idx in enumerate(OrderedSet(sorted_asu_indices)):
            n_in_group = self.h_index_matrix.col(j).non_zeroes
            if miller_idx in target_asu_Ih_dict:
                i = location_in_unscaled_array
                new_Ih_values.set_selected(
                    flex.size_t(range(i, i + n_in_group)),
                    flex.double(n_in_group, target_asu_Ih_dict[miller_idx]),
                )
            location_in_unscaled_array += n_in_group
        self.Ih_values.set_selected(permuted, new_Ih_values)
        sel = self.Ih_values != 0.0
        new_table = self.select(sel)
        # now set attributes to update object
        self.Ih_table = new_table.Ih_table
        self.h_index_matrix = new_table.h_index_matrix
        self.h_expand_matrix = new_table.h_expand_matrix
        self.block_selections = new_table.block_selections

    @property
    def inverse_scale_factors(self):
        """The inverse scale factors of the reflections."""
        return self.Ih_table["inverse_scale_factor"]

    @inverse_scale_factors.setter
    def inverse_scale_factors(self, new_scales):
        if new_scales.size() != self.size:
            assert 0, """attempting to set a new set of scale factors of different
      length than previous assignment: was %s, attempting %s""" % (
                self.inverse_scale_factors.size(),
                new_scales.size(),
            )
        else:
            self.Ih_table["inverse_scale_factor"] = new_scales

    @property
    def variances(self):
        """The variances of the reflections."""
        return self.Ih_table["variance"]

    @property
    def intensities(self):
        """The unscaled reflection intensities."""
        return self.Ih_table["intensity"]

    @property
    def Ih_values(self):
        """The bset-estimated intensities of symmetry equivalent reflections."""
        return self.Ih_table["Ih_values"]

    @property
    def weights(self):
        """The weights that will be used in scaling."""
        return self.Ih_table["weights"]

    @weights.setter
    def weights(self, new_weights):
        if new_weights.size() != self.size:
            assert 0, """attempting to set a new set of weights of different
      length than previous assignment: was %s, attempting %s""" % (
                self.size,
                new_weights.size(),
            )
        self.Ih_table["weights"] = new_weights

    @property
    def size(self):
        """Return the length of the stored Ih_table (a reflection table)."""
        return self.Ih_table.size()

    @property
    def asu_miller_index(self):
        """Return the miller indices in the asymmetric unit."""
        return self.Ih_table["asu_miller_index"]

    def setup_binner(self, unit_cell, space_group, n_resolution_bins):
        ma = _reflection_table_to_iobs(self.Ih_table, unit_cell, space_group)
        # need d star sq step
        d_star_sq = ma.d_star_sq().data()
        d_star_sq_min = flex.min(d_star_sq)
        d_star_sq_max = flex.max(d_star_sq)
        span = d_star_sq_max - d_star_sq_min
        relative_tolerance = 1e-6
        d_star_sq_max += span * relative_tolerance
        d_star_sq_min -= span * relative_tolerance
        step = (d_star_sq_max - d_star_sq_min) / n_resolution_bins

        self.binner = ma.setup_binner_d_star_sq_step(
            auto_binning=False,
            d_max=uctbx.d_star_sq_as_d(d_star_sq_max),
            d_min=uctbx.d_star_sq_as_d(d_star_sq_min),
            d_star_sq_step=step,
        )


def _reflection_table_to_iobs(table, unit_cell, space_group):

    miller_set = miller.set(
        crystal_symmetry=crystal.symmetry(
            unit_cell=unit_cell,
            space_group=space_group,
            assert_is_compatible_unit_cell=False,
        ),
        indices=table["asu_miller_index"],
        anomalous_flag=False,
    )
    i_obs = miller.array(
        miller_set, data=table["intensity"] / table["inverse_scale_factor"]
    )
    i_obs.set_observation_type_xray_intensity()
    i_obs.set_sigmas((table["variance"] ** 0.5) / table["inverse_scale_factor"])
    i_obs.set_info(miller.array_info(source="DIALS", source_type="reflection_tables"))
    return i_obs

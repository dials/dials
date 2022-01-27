"""
Algorithm to select a well connected subset of reflections for scaling.

Description of reflection selection algorithm. To get the best 'connectedness',
we want to select groups of reflections which belong to more than one class.
A class in this case is a volume of reciprocal space (e.g 12 areas on the
surface of a sphere) or a class can be a dataset (i.e. n classes for n datasets).

First we construct the following matrix of classes vs symmetry groups. For
example, this matrix describes 3 datasets with 7 symmetry unique groups.

              symmetry groups
           0  1  2  3  4  5  6
        0  3  3  2  0  1  1  1
classes 1  0  2  0  0  3  2  1
        2  2  1  1  5  0  4  0

Here the entries of the matrix are the number of reflections belonging to the
group and class. Then, the matrix is sorted by the number of classes that each
group covers, i.e. the number of nonzero entries in the column:
number of nonzero entries: [2, 3, 2, 1, 2, 3, 2]
sorted matrix:
              symmetry groups
           1  5  0  2  4  6  3
        0  3  1  3  2  1  1  0
classes 1  2  2  0  0  3  1  0
        2  1  4  2  1  0  0  5

Now, we choose a target number of reflections per class e.g. 5. To choose the
reflection groups, we start with the first column.
number of chosen reflections per class: [3, 2, 1]
symmetry groups used:                   [1]

To determine the next group to add, we search for the first group (matrix column)
that has a reflection in the least populated class so far i.e. class 2.
In this case, the first unused group is group 5:
number of chosen reflections per class: [4, 4, 5]
symmetry groups used:                   [1, 5]
In this way, we build up the dataset by choosing the highest-connected groups
that have a reflection in the most-deficient class.

Next we need to add a group with a reflection in class 0 (first we find is group 0):
number of chosen reflections per class: [7, 4, 7]
symmetry groups used:                   [1, 5, 0]

Next we need to add a group with a reflection in class 1 (first we find is group 4):
number of chosen reflections per class: [8, 7, 7]
symmetry groups used:                   [1, 5, 0, 4]

We have now reached our target for all classes and can therefore stop.
The four symmetry groups are the highest connected groups that give use good
coverage across all classes, and are therefore the best reflections to use for
minimisation. If there were fewer reflections in one class than the target,
then this algorithm will add all groups with reflections in that class and then
continue with the remaining classes.

For single dataset minimisation, this algorithm is used to select reflection
groups with good reciprocal space coverage, repeated across resolution bins.
For multi dataset minimisation, this algorithm is also used to select highly
connected reflections between datasets. The reflections used for minimisation
are those which are selected by either method - inter-dataset connectedness or
intra-dataset connectedness.
"""

from __future__ import annotations

import logging
from math import floor

import numpy as np

from dxtbx import flumpy
from scitbx import sparse

from dials.algorithms.scaling.scaling_utilities import (
    BadDatasetForScalingException,
    Reasons,
)
from dials.array_family import flex
from dials.util import tabulate

logger = logging.getLogger("dials")


def _build_class_matrix(class_index, class_matrix, offset=0):
    for (i, val) in enumerate(class_index, start=offset):
        class_matrix[val, i] = 1.0
    return class_matrix


def _select_groups_on_Isigma_cutoff(Ih_table, cutoff=2.0):
    """Select groups with multiplicity>1, Isigma>cutoff"""
    sumIsigm = Ih_table.sum_in_groups(
        Ih_table.intensities / np.sqrt(Ih_table.variances)
    )
    n = Ih_table.group_multiplicities()
    avg_Isigma = sumIsigm / n
    sel = avg_Isigma > cutoff
    sel2 = n > 1
    if not sel2.any():
        raise SystemExit(
            """
Could not find any cross-dataset connected reflections with multiplicity > 1,
scaling not possible."""
        )
    sel &= sel2
    if not sel.any():
        logger.warning(
            """
Warning: Could not select any reflections for <I/sI> > %s.
Reducing Isigma_cutoff to zero to attempt continuation.""",
            cutoff,
        )
        sel = avg_Isigma > 0.0 & sel2
        if not sel.any():
            raise SystemExit(
                """
Could not find any cross-dataset connected groups with <I/sI> > 0,
scaling not possible."""
            )
    sel_Ih_table = Ih_table.select_on_groups(sel)
    return sel_Ih_table


def _perform_quasi_random_selection(
    Ih_table, n_datasets, min_per_class, min_total, max_total
):

    class_matrix = sparse.matrix(n_datasets, Ih_table.size)
    class_matrix = _build_class_matrix(
        flumpy.from_numpy(Ih_table.Ih_table["dataset_id"].to_numpy()), class_matrix
    )
    segments_in_groups = class_matrix * Ih_table.h_index_matrix
    total = flex.double(segments_in_groups.n_cols, 0)
    for i, col in enumerate(segments_in_groups.cols()):
        total[i] = col.non_zeroes
    perm = flex.sort_permutation(total, reverse=True)
    sorted_class_matrix = segments_in_groups.select_columns(perm)
    # matrix of segment index vs asu groups

    # now want to fill up until good coverage across board
    total_in_classes, cols_not_used = _loop_over_class_matrix(
        sorted_class_matrix, min_per_class, min_total, max_total
    )

    cols_used = flex.bool(sorted_class_matrix.n_cols, True)
    cols_used.set_selected(cols_not_used, False)
    actual_cols_used = perm.select(cols_used)

    # now need to get reflection selection
    reduced_Ih = Ih_table.select_on_groups(actual_cols_used)
    indices_this_res = reduced_Ih.Ih_table["loc_indices"]
    dataset_ids_this_res = reduced_Ih.Ih_table["dataset_id"]

    n_groups_used = len(actual_cols_used)

    return (
        flumpy.from_numpy(indices_this_res.to_numpy()),
        flumpy.from_numpy(dataset_ids_this_res.to_numpy()),
        n_groups_used,
        total_in_classes,
    )


def select_connected_reflections_across_datasets(
    Ih_table, experiment, Isigma_cutoff=2.0, min_total=40000, n_resolution_bins=20
):
    """Select highly connected reflections across datasets."""
    assert Ih_table.n_work_blocks == 1
    Ih_table = Ih_table.Ih_table_blocks[0]
    sel_Ih_table = _select_groups_on_Isigma_cutoff(Ih_table, Isigma_cutoff)

    # now split into resolution bins
    sel_Ih_table.setup_binner(
        experiment.crystal.get_unit_cell(),
        experiment.crystal.get_space_group(),
        n_resolution_bins,
    )
    binner = sel_Ih_table.binner

    # prepare parameters for selection algorithm.
    n_datasets = len(set(sel_Ih_table.Ih_table["dataset_id"].to_numpy()))
    min_per_class = min_total / (n_datasets * 4.0)
    max_total = min_total * 1.2
    logger.info(
        """
Using quasi-random reflection selection. Selecting from %s symmetry groups
with <I/sI> > %s (%s reflections)). Selection target of %.2f reflections
from each dataset, with a total number between %.2f and %.2f.
""",
        sel_Ih_table.n_groups,
        Isigma_cutoff,
        sel_Ih_table.size,
        min_per_class,
        min_total,
        max_total,
    )
    # split across resolution bins
    mpc = int(min_per_class / n_resolution_bins)
    mint = int(min_total / n_resolution_bins)
    maxt = int(max_total / n_resolution_bins)

    header = ["d-range", "n_groups", "n_refl"] + [str(i) for i in range(n_datasets)]
    rows = []
    if n_datasets >= 15:
        summary_rows = []
        summary_header = ["d-range", "n_groups", "n_refl"]

    indices = flex.size_t()
    dataset_ids = flex.size_t()
    total_groups_used = 0
    n_cols_used = 0

    for ibin in binner.range_all():
        sel = binner.selection(ibin)
        res_Ih_table = sel_Ih_table.select(flumpy.to_numpy(sel))
        if not res_Ih_table.Ih_table.size:
            continue

        (
            indices_this_res,
            dataset_ids_this_res,
            n_groups_used,
            total_per_dataset,
        ) = _perform_quasi_random_selection(res_Ih_table, n_datasets, mpc, mint, maxt)

        indices.extend(indices_this_res)
        dataset_ids.extend(dataset_ids_this_res)
        total_groups_used += n_groups_used
        d0, d1 = binner.bin_d_range(ibin)
        drange = str(round(d0, 3)) + " - " + str(round(d1, 3))
        n_refl = str(int(indices_this_res.size()))
        rows.append(
            [drange, str(n_groups_used), n_refl]
            + [str(int(i)) for i in total_per_dataset]
        )
        if n_datasets >= 15:
            summary_rows.append([drange, str(n_groups_used), n_refl])
        n_cols_used += n_groups_used

    logger.info(
        "Summary of cross-dataset reflection groups chosen (%s groups, %s reflections):",
        n_cols_used,
        indices.size(),
    )
    if n_datasets < 15:
        logger.info(tabulate(rows, header))
    else:
        logger.info(tabulate(summary_rows, summary_header))
        logger.debug(tabulate(rows, header))

    return indices, dataset_ids


def _loop_over_class_matrix(
    sorted_class_matrix, min_per_area, min_per_bin, max_per_bin
):
    """Build up the reflection set by looping over the class matrix."""

    def _get_next_row_needed(total_in_classes):
        current_min = flex.min(total_in_classes)
        for i, val in enumerate(total_in_classes):
            if val == current_min:
                row_needed = i
                break
        return row_needed

    def _add_next_column(
        cols_not_used, row_needed, sorted_class_matrix, total_in_classes
    ):
        for i, col in enumerate(cols_not_used):
            if sorted_class_matrix.col(col)[row_needed] != 0.0:
                total_in_classes += sorted_class_matrix.col(col).as_dense_vector()
                del cols_not_used[i]
                return cols_not_used, total_in_classes, True
        # else couldn't find enough of this one!
        return cols_not_used, total_in_classes, False

    total_in_classes = sorted_class_matrix.col(0).as_dense_vector()
    defecit = flex.double(sorted_class_matrix.n_rows, 0)
    cols_not_used = flex.size_t(range(1, sorted_class_matrix.n_cols))
    total_deficit = 0
    while (
        flex.min(total_in_classes) < min_per_area
        and (flex.sum(total_in_classes) - total_deficit) < max_per_bin
    ):
        # first find which class need most of
        row_needed = _get_next_row_needed(total_in_classes)
        # now try to add the most-connected column that includes that class
        cols_not_used, total_in_classes, success = _add_next_column(
            cols_not_used, row_needed, sorted_class_matrix, total_in_classes
        )
        # return whether successful, updated totals and which cols are left.
        if not success:
            # want to stop looking for that class as no more left
            current_in_row = total_in_classes[row_needed]
            defecit[row_needed] = min_per_area - current_in_row
            total_deficit += min_per_area - current_in_row
            total_in_classes[row_needed] = min_per_area
        if flex.sum(total_in_classes) > max_per_bin:
            # if we have reached the maximum, then finish there
            return total_in_classes - defecit, cols_not_used
    total_in_classes -= defecit
    n = flex.sum(total_in_classes)
    # if we haven't reached the minimum total, then need to add more until we
    # reach it or run out of reflections
    if n < min_per_bin and cols_not_used:
        # how many have deficit? (i.e. no more left?)
        c = sum(1 for d in defecit if d != 0.0)
        n_classes = sorted_class_matrix.n_rows
        multiplier = int(floor(min_per_bin * (n_classes - c) / (n * n_classes)) + 1)
        new_limit = min_per_area * multiplier  # new limit per area

        for i, d in enumerate(defecit):
            if d != 0.0:
                # don't want to be searching for those classes that we know dont have any left
                total_in_classes[i] = new_limit
                defecit[i] = d + new_limit - min_per_area
        while cols_not_used and flex.min(total_in_classes) < new_limit:
            row_needed = _get_next_row_needed(total_in_classes)
            cols_not_used, total_in_classes, success = _add_next_column(
                cols_not_used, row_needed, sorted_class_matrix, total_in_classes
            )
            if not success:
                current_in_row = total_in_classes[row_needed]
                defecit[row_needed] = new_limit - current_in_row
                total_in_classes[row_needed] = new_limit
        return total_in_classes - defecit, cols_not_used
    return total_in_classes, cols_not_used


def _determine_Isigma_selection(reflection_table, params):
    Ioversigma = reflection_table["intensity"] / flex.sqrt(reflection_table["variance"])
    Isiglow, Isighigh = params.reflection_selection.Isigma_range
    selection = Ioversigma > Isiglow
    if Isighigh != 0.0:
        selection &= Ioversigma < Isighigh
        reason = f"in I/sigma range ({Isighigh} > I/sig > {Isiglow})"
    else:
        reason = f"in I/sigma range (I/sig > {Isiglow})"
    return selection, reason


def _determine_partiality_selection(reflection_table, params):
    min_partiality = params.reflection_selection.min_partiality
    selection = reflection_table["partiality"] > min_partiality
    reason = f"above min partiality ( > {min_partiality})"
    return selection, reason


def _determine_d_range_selection(reflection_table, params):
    d_min, d_max = params.reflection_selection.d_range
    d_sel = reflection_table["d"] > d_min
    d_sel &= reflection_table["d"] < d_max
    reason = f"in d range ({d_max} > d > {d_min})"
    return d_sel, reason


def _determine_E2_range_selection(reflection_table, params):
    Elow, Ehigh = params.reflection_selection.E2_range
    sel1 = reflection_table["Esq"] > Elow
    sel2 = reflection_table["Esq"] < Ehigh
    Esq_sel = sel1 & sel2
    reason = f"in E^2 range ({Ehigh} > E^2 > {Elow})"
    return Esq_sel, reason


def calculate_scaling_subset_ranges(reflection_table, params, print_summary=False):
    selection, reasons = _common_range_selections(Reasons(), reflection_table, params)
    if print_summary:
        logger.info(
            "%s reflections were preselected for scale factor determination \n"
            + "out of %s suitable reflections: \n%s",
            selection.count(True),
            reflection_table.size(),
            reasons,
        )
    if selection.count(True) == 0:
        raise BadDatasetForScalingException(
            """No reflections pass all user-controllable selection criteria"""
        )
    return selection


def _common_range_selections(reasons, reflection_table, params):
    selection, reason = _determine_Isigma_selection(reflection_table, params)
    reasons.add_reason(reason, selection.count(True))
    if "partiality" in reflection_table:
        sel, reason = _determine_partiality_selection(reflection_table, params)
        reasons.add_reason(reason, sel.count(True))
        selection &= sel
    if params.reflection_selection.d_range:
        sel, reason = _determine_d_range_selection(reflection_table, params)
        reasons.add_reason(reason, sel.count(True))
        selection &= sel
    return selection, reasons


def calculate_scaling_subset_ranges_with_E2(reflection_table, params):
    """Select reflections with non-zero weight and update scale weights."""
    reasons = Reasons()
    selection = ~reflection_table.get_flags(
        reflection_table.flags.user_excluded_in_scaling
    )
    selection &= ~reflection_table.get_flags(
        reflection_table.flags.excluded_for_scaling
    )
    reasons.add_reason("suitable/selected for scaling", selection.count(True))
    if reflection_table["Esq"].count(1.0) != reflection_table.size():
        sel, reason = _determine_E2_range_selection(reflection_table, params)
        reasons.add_reason(reason, sel.count(True))
        selection &= sel
    sel, reasons = _common_range_selections(reasons, reflection_table, params)
    selection &= sel
    logger.info(
        "%s reflections were selected for scale factor determination \n"
        + "out of %s suitable reflections: \n%s",
        selection.count(True),
        reflection_table.size(),
        reasons,
    )
    if selection.count(True) == 0:
        raise BadDatasetForScalingException(
            """No reflections pass all user-controllable selection criteria"""
        )
    return selection

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

Now, we chose a target number of reflections per class e.g. 5. To chose the
reflection groups, we start with the first column.
number of chosen reflections per class: [3, 2, 1]
symmetry groups used:                   [1]

To determine the next group to add, we search for the first group (matrix column)
that has a reflection in the least populated class so far i.e. class 2.
In this case, the first unused group is group 5:
number of chosen reflections per class: [4, 4, 5]
symmetry groups used:                   [1, 5]
In this way, we build up the dataset by chosing the highest-connected groups
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
from __future__ import absolute_import, division, print_function
import logging
from math import pi, floor
import libtbx
from libtbx.table_utils import simple_table
from dials.array_family import flex
from dials.algorithms.scaling.scaling_utilities import (
    Reasons,
    BadDatasetForScalingException,
)
from dials_scaling_ext import calc_theta_phi

logger = logging.getLogger("dials")

magic_theta = 70.53  # arccos(1/3) - to make all sphere regions equal.


def determine_reflection_selection_parameters(params, experiments, reflections):
    """Algorithm to determine suitable parameters for the quasi-random algorithm.

    Only triggered when reflection_selection.method == auto
    Rules:
      if < 20000 total reflections, just use them all.
      for single dataset minimisation:
        - want at least the greater of 100 * n_params / 2% of reflections, as long
          as this is between 10000 and 50000 reflections.
        - use 20 resolution bins
      for multi-dataset minimisation:
        - try to get at least 100-200 reflections per parameter
        - try to get a split of 20% inter-dataset, 80% intra-dataset reflections
        - try to use at least 50000 reflections (no upper limit)
        - if trying to use <= 360 reflections per dataset, use 10 resolution bins,
          else use 20 resolution bins (to get at least 2 refl per area).
    """
    if params.reflection_selection.method in (None, libtbx.Auto, "auto"):
        if len(experiments) == 1:
            n_suitable_refl = (
                reflections[0]
                .get_flags(reflections[0].flags.bad_for_scaling, all=False)
                .count(False)
            )
            if n_suitable_refl < 20000:
                params.reflection_selection.method = "use_all"
                logger.info(
                    "Using all reflections for minimisation as less than 20000 suitable reflections."
                )
                return
            n_params = experiments[0].scaling_model.n_params
            num = 100 * n_params  # want at least this
            two_percent = 0.02 * n_suitable_refl
            n_to_use = max(10000, num, two_percent)  # in reality,
            # allows up to ~100K due to overfilling effect?
            # but with 12 areas x 20 res bins = 240
            n_resolution_bins = 20
            min_per_area = int(n_to_use / 240.0)  # (42 to 412 per area)
            params.reflection_selection.method = "quasi_random"
            params.reflection_selection.quasi_random.min_per_area = [min_per_area]
            params.reflection_selection.quasi_random.n_resolution_bins = [
                n_resolution_bins
            ]
            logger.info(
                """Determined quasi-random reflection selection parameters:
min_per_area : %s, n_resolution_bins: %s""",
                min_per_area,
                n_resolution_bins,
            )
        else:
            if sum([r.size() for r in reflections]) < 40000:
                params.reflection_selection.method = "use_all"
                logger.info(
                    "Using all reflections for minimisation as less than 40000 reflections in total."
                )
                return
            # this attempts to get ~ 100-200 refl per param
            n_params = [exp.scaling_model.n_params for exp in experiments]
            min_per_class = 120 * sum(n_params) / len(n_params)
            total_target = 2.0 * 80 * sum(n_params) + (min_per_class * len(n_params))
            # extra factor of 2 as set min, max limits of 1.5 to 3.0 - so get 40/160 split
            min_per_area = []
            n_resolution_bins = []
            scale = 1.0
            if total_target < 50000:
                scale = 50000 / total_target
                min_per_class = int(min_per_class * scale)
            else:
                min_per_class = int(min_per_class)
            for n_param in n_params:
                target = 80 * n_param * scale
                if target <= 360:
                    n_resolution_bins.append(10)
                    min_per_area.append(int(target / 120.0))
                else:
                    n_resolution_bins.append(20)
                    min_per_area.append(int(target / 240.0))
            params.reflection_selection.method = "quasi_random"
            params.reflection_selection.quasi_random.min_per_area = min_per_area
            params.reflection_selection.quasi_random.n_resolution_bins = (
                n_resolution_bins
            )
            params.reflection_selection.quasi_random.multi_dataset.min_per_dataset = (
                min_per_class
            )
            logger.debug(
                """Determined quasi-random reflection selection parameters:
min_per_area: \n%s\nn_resolution_bins: \n%s\nmin_per_class %s""",
                min_per_area,
                n_resolution_bins,
                min_per_class,
            )
    elif params.reflection_selection.method == "quasi_random" and len(reflections) > 1:
        mpa = params.reflection_selection.quasi_random.min_per_area
        nrb = params.reflection_selection.quasi_random.n_resolution_bins
        if len(mpa) == 1:
            mpa *= len(reflections)
        elif len(mpa) != len(reflections):
            inital = mpa[0]
            params.reflection_selection.quasi_random.min_per_area = [inital] * len(
                reflections
            )
            logger.warning(
                """Warning:
Using quasi-random reflection selection with manual parameters, but length
of min_per_area list (%s) not equal to number of reflection tables (%s).
Using first min_per_area value for all datasets.\n""",
                len(mpa),
                len(reflections),
            )
        if len(nrb) == 1:
            nrb *= len(reflections)
        elif len(nrb) != len(reflections):
            initial = nrb[0]
            params.reflection_selection.quasi_random.n_resolution_bins = [
                initial
            ] * len(reflections)
            logger.warning(
                """Warning:
Using quasi-random reflection selection with manual parameters, but length
of n_resolution_bins list (%s) not equal to number of reflection tables (%s).
Using first n_resolution_bins value for all datasets.\n""",
                len(nrb),
                len(reflections),
            )


def assign_segment_index(reflections):
    """Divide a sphere into 12 segments of equal area."""
    seg_idx = flex.int(reflections.size(), -1)
    all_theta = flex.bool(reflections.size(), True)
    theta_idx_0 = reflections["theta"] < magic_theta
    theta_idx_2 = reflections["theta"] > (180.0 - magic_theta)
    theta_idx_1 = all_theta & ~(theta_idx_0 | theta_idx_2)

    phi_norm = flex.floor(reflections["phi"] / 90.0)
    phi_iround = phi_norm.iround()
    phi_idx = phi_iround % 4  # phi from -180 to 180
    seg_idx.set_selected(theta_idx_0, phi_idx)
    phi_idx += 8
    seg_idx.set_selected(theta_idx_2, phi_idx)
    phi_idx = (flex.floor((reflections["phi"] + 45.0) / 90.00001).iround() % 4) + 4
    seg_idx.set_selected(theta_idx_1, phi_idx)
    assert flex.min(seg_idx) >= 0, flex.min(seg_idx)
    assert flex.max(seg_idx) <= 11, flex.max(seg_idx)
    reflections["class_index"] = seg_idx
    return reflections


def _build_class_matrix(reflections, class_matrix, offset=0):
    for (i, val) in enumerate(reflections["class_index"], start=offset):
        class_matrix[val, i] = 1.0
    return class_matrix


def select_highly_connected_reflections(
    Ih_table_block, experiment, min_per_area, n_resolution_bins, print_summary=False
):
    """Select highly connected reflections within a dataset, across resolutions."""
    min_per_bin = min_per_area * 12 * 1.5
    max_per_bin = min_per_area * 12 * 3.0
    assert "s1c" in Ih_table_block.Ih_table

    theta_phi_1 = calc_theta_phi(Ih_table_block.Ih_table["s1c"])
    theta = theta_phi_1.parts()[0]
    phi = theta_phi_1.parts()[1]
    Ih_table_block.Ih_table["phi"] = (phi * 180 / pi) + 180.0
    Ih_table_block.Ih_table["theta"] = theta * 180 / pi

    Ih_table_block.Ih_table = assign_segment_index(Ih_table_block.Ih_table)
    Ih_table_block.setup_binner(
        experiment.crystal.get_unit_cell(),
        experiment.crystal.get_space_group(),
        n_resolution_bins,
    )
    binner = Ih_table_block.binner

    overall_indices = flex.size_t()

    header = ["d-range", "n_refl"] + [str(i) for i in range(0, 12)]
    rows = []

    for ibin in binner.range_all():
        sel = binner.selection(ibin)
        sel_Ih_table_block = Ih_table_block.select(sel)
        indices_wrt_original = Ih_table_block.Ih_table["loc_indices"].select(sel)
        indices, total_in_classes = select_highly_connected_reflections_in_bin(
            sel_Ih_table_block, min_per_area, min_per_bin, max_per_bin
        )
        if indices:
            overall_indices.extend(indices_wrt_original.select(indices))
            d0, d1 = binner.bin_d_range(ibin)
            rows.append(
                [
                    str(round(d0, 3)) + " - " + str(round(d1, 3)),
                    str(int(flex.sum(total_in_classes))),
                ]
                + [str(int(i)) for i in total_in_classes]
            )
    st = simple_table(rows, header)
    msg = """\nSummary of reflection selection algorithm for this dataset:
%s resolution bins, target: at least %s reflections per area,
between %s and %s reflections per resolution bin""" % (
        n_resolution_bins,
        min_per_area,
        18 * min_per_area,
        36 * min_per_area,
    )
    if print_summary:
        logger.info(msg)
        logger.info(st.format())
    else:
        logger.debug(msg)
        logger.debug(st.format())
    return overall_indices


def select_connected_reflections_across_datasets(
    Ih_table, min_per_class=500, min_multiplicity=2, Isigma_cutoff=2.0
):
    """Select highly connected reflections across datasets."""
    assert Ih_table.n_work_blocks == 1
    Ih_table = Ih_table.Ih_table_blocks[0]
    n_refl = Ih_table.Ih_table.size()
    # first select groups where avg I/sigma > cutoff?
    I_over_sigma = Ih_table.intensities / (Ih_table.variances ** 0.5)
    sumIsigm = I_over_sigma * Ih_table.h_index_matrix
    n = flex.double(I_over_sigma.size(), 1.0) * Ih_table.h_index_matrix
    avg_Isigma = sumIsigm / n
    sel = avg_Isigma > Isigma_cutoff
    sel2 = n >= min_multiplicity
    sel &= sel2
    if sel.count(True) == 0:
        if sel2.count(True) == 0:
            raise SystemExit(
                "Could not find any cross-dataset connected reflections with min_multiplicity >= %s."
                % min_multiplicity
            )
        logger.warning(
            """Warning:
Could not select any reflections for <I/sI> > %s and min_multiplicity >= %s.
Reducing Isigma_cutoff to zero to attempt continuation.""",
            Isigma_cutoff,
            min_multiplicity,
        )
        sel = avg_Isigma > 0.0
        sel &= sel2
        if sel.count(True) == 0:
            raise SystemExit(
                """Could not find any cross-dataset connected reflections belonging
to groups with for <I/sI> > 0 for min_multiplicity >= %s"""
                % min_multiplicity
            )

    sel_Ih_table = Ih_table.select_on_groups(sel)
    logger.info(
        """
Determining highly connected reflections across datasets for scaling model
minimisation. Prefiltering for symmetry groups that have reflections with
multiplicity of at least %s with <I/sI> > %s. This leaves %s out of %s reflections
to potentially use for scaling model minimisation based on cross-dataset
connectedness (these belong to %s symmetry groups).""",
        min_multiplicity,
        Isigma_cutoff,
        sel_Ih_table.Ih_table.size(),
        n_refl,
        sel.count(True),
    )

    n_datasets = len(set(sel_Ih_table.Ih_table["dataset_id"]))
    min_total = min_per_class * n_datasets
    max_total = min_total * 3.0

    logger.info(
        """
Attempting to choose at least %s reflections from each dataset,
with a total number between %s and %s.""",
        min_per_class,
        min_total,
        max_total,
    )

    from scitbx import sparse

    class_matrix = sparse.matrix(n_datasets, sel_Ih_table.Ih_table.size())

    sel_Ih_table.Ih_table["class_index"] = sel_Ih_table.Ih_table["dataset_id"]

    class_matrix = _build_class_matrix(sel_Ih_table.Ih_table, class_matrix)
    segments_in_groups = class_matrix * sel_Ih_table.h_index_matrix

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
    reduced_Ih = sel_Ih_table.select_on_groups_isel(actual_cols_used)
    indices = reduced_Ih.Ih_table["loc_indices"]
    dataset_ids = reduced_Ih.Ih_table["dataset_id"]
    logger.info(
        """
Choosing %s cross-dataset connected reflections from %s symmetry groups for minimisation.\n""",
        indices.size(),
        len(actual_cols_used),
    )
    return indices, dataset_ids, total_in_classes


def select_highly_connected_reflections_in_bin(
    Ih_table_block, min_per_class=2, min_total=1000, max_total=10000
):
    """Select highly connected reflections within a resolution shell."""

    n = flex.double(Ih_table_block.size, 1.0) * Ih_table_block.h_index_matrix
    sel = n > 1
    if sel.count(True) == 0:
        return None, None

    Ih_table_block.Ih_table["loc_indices"] = flex.size_t(range(0, Ih_table_block.size))
    Ih_table_block = Ih_table_block.select_on_groups(sel)

    from scitbx import sparse

    class_matrix = sparse.matrix(12, Ih_table_block.size)

    class_matrix = _build_class_matrix(Ih_table_block.Ih_table, class_matrix)
    segments_in_groups = class_matrix * Ih_table_block.h_index_matrix

    total = flex.int(segments_in_groups.n_cols, 0)
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
    reduced_Ih = Ih_table_block.select_on_groups_isel(actual_cols_used)
    indices = reduced_Ih.Ih_table["loc_indices"]
    return indices, total_in_classes


def _get_next_row_needed(total_in_classes):
    current_min = flex.min(total_in_classes)
    for i, val in enumerate(total_in_classes):
        if val == current_min:
            row_needed = i
            break
    return row_needed


def _add_next_column(cols_not_used, row_needed, sorted_class_matrix, total_in_classes):
    for i, col in enumerate(cols_not_used):
        if sorted_class_matrix.col(col)[row_needed] != 0.0:
            total_in_classes += sorted_class_matrix.col(col).as_dense_vector()
            del cols_not_used[i]
            return cols_not_used, total_in_classes, True
    # else couldn't find enough of this one!
    return cols_not_used, total_in_classes, False


def _loop_over_class_matrix(
    sorted_class_matrix, min_per_area, min_per_bin, max_per_bin
):
    """Build up the reflectio set by looping over the class matrix."""
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
        multiplier = int(floor(min_per_bin / n) + 1)
        new_limit = min_per_area * multiplier
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


def calculate_scaling_subset_connected(
    global_Ih_table, experiment, params, preselection=None, print_summary=False
):
    """Determine the selection for the reflection table of suitable reflections.

    A preselection can be given, which is applied to the global_Ih_table before
    determining the connected subset."""
    min_per_area = params.reflection_selection.quasi_random.min_per_area[0]
    n_resolution_bins = params.reflection_selection.quasi_random.n_resolution_bins[0]
    suitable_selection = flex.bool(global_Ih_table.size, False)
    if preselection:
        Ih_table = global_Ih_table.select(preselection)
    else:
        Ih_table = global_Ih_table
    indices = select_highly_connected_reflections(
        Ih_table, experiment, min_per_area, n_resolution_bins, print_summary
    )
    suitable_selection.set_selected(indices, True)
    if print_summary:
        logger.info(
            "%s reflections were selected for scale factor determination \n"
            + "out of %s suitable reflections. ",
            suitable_selection.count(True),
            global_Ih_table.size,
        )
    if suitable_selection.count(True) == 0:
        raise BadDatasetForScalingException(
            """No reflections pass all user-controllable selection criteria"""
        )
    return suitable_selection


def _determine_Isigma_selection(reflection_table, params):
    Ioversigma = reflection_table["intensity"] / (reflection_table["variance"] ** 0.5)
    Isiglow, Isighigh = params.reflection_selection.Isigma_range
    selection = Ioversigma > Isiglow
    if Isighigh != 0.0:
        selection &= Ioversigma < Isighigh
        reason = "in I/sigma range (%s > I/sig > %s)" % (Isighigh, Isiglow)
    else:
        reason = "in I/sigma range (I/sig > %s)" % Isiglow
    return selection, reason


def _determine_partiality_selection(reflection_table, params):
    min_partiality = params.reflection_selection.min_partiality
    selection = reflection_table["partiality"] > min_partiality
    reason = "above min partiality ( > %s)" % min_partiality
    return selection, reason


def _determine_d_range_selection(reflection_table, params):
    d_min, d_max = params.reflection_selection.d_range
    d_sel = reflection_table["d"] > d_min
    d_sel &= reflection_table["d"] < d_max
    reason = "in d range (%s > d > %s)" % (d_max, d_min)
    return d_sel, reason


def _determine_E2_range_selection(reflection_table, params):
    Elow, Ehigh = params.reflection_selection.E2_range
    sel1 = reflection_table["Esq"] > Elow
    sel2 = reflection_table["Esq"] < Ehigh
    Esq_sel = sel1 & sel2
    reason = "in E^2 range (%s > E^2 > %s)" % (Ehigh, Elow)
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

"""Calculations relevant to reporting."""
from __future__ import absolute_import, division, print_function

from cctbx import miller
from dials.algorithms.scaling.scaling_library import scaled_data_as_miller_array
from dials.util.batch_handling import (
    calculate_batch_offsets,
    get_batch_ranges,
    assign_batches_to_reflections,
)
from libtbx.str_utils import make_sub_header
from scitbx.array_family import flex
from six.moves import cStringIO as StringIO


def batch_dependent_properties(batches, intensities, scales=None):
    """
    Calculate intensity properties as a function of batch.

    Manage individual calculations to ensure that all calculations are
    done on the same batches array and consistent length lists returned.

    Args:
        batches (miller array): the batch numbers for the reflections.
        intensities (miller array): the reflection intensities, with
            sigmas also present.
        scales (miller array, optional): the scale factors of the reflections.

    Returns:
        (tuple): tuple containing:
            binned_batches (list): list of batch number bins
            rmerge (list): rmerge for each bin
            isigi (list): average I/sigma for each bin
            scales (list): average scale for each bin if scales given, else None

    Raises:
        AssertionError: Raised if sigmas not present in intensities, or if
            arrays sizes do not match.

    """
    assert intensities.size() == batches.size()
    if scales:
        assert scales.size() == batches.size()
    assert intensities.sigmas() is not None
    sel = intensities.sigmas() > 0
    intensities = intensities.select(sel)
    batches = batches.select(sel)
    if scales:
        scales = scales.select(sel)

    binned_batches, rmerge = rmerge_vs_batch(intensities, batches)
    _, isigi = i_sig_i_vs_batch(intensities, batches)

    scalesvsbatch = None
    if scales:
        _, scalesvsbatch = scales_vs_batch(scales, batches)
    return binned_batches, rmerge, isigi, scalesvsbatch


def combined_table_to_batch_dependent_properties(
    combined_table, experiments, scaled_array=None
):
    """Extract batch dependent properties from a combined reflection table."""
    tables = []
    for id_ in set(combined_table["id"]).difference({-1}):
        tables.append(combined_table.select(combined_table["id"] == id_))

    return reflection_tables_to_batch_dependent_properties(
        tables, experiments, scaled_array
    )


def reflection_tables_to_batch_dependent_properties(
    reflection_tables, experiments, scaled_array=None
):
    """Extract batch dependent properties from a reflection table list."""
    offsets = calculate_batch_offsets(experiments)
    reflection_tables = assign_batches_to_reflections(reflection_tables, offsets)
    # filter bad refls and negative scales
    batches = flex.int()
    scales = flex.double()
    for r in reflection_tables:
        sel = ~r.get_flags(r.flags.bad_for_scaling, all=False)
        sel &= r["inverse_scale_factor"] > 0
        batches.extend(r["batch"].select(sel))
        scales.extend(r["inverse_scale_factor"].select(sel))
    if not scaled_array:
        scaled_array = scaled_data_as_miller_array(reflection_tables, experiments)
    ms = scaled_array.customized_copy()
    batch_array = miller.array(ms, data=batches)

    batch_ranges = get_batch_ranges(experiments, offsets)
    batch_data = [{"id": i, "range": r} for i, r in enumerate(batch_ranges)]

    properties = batch_dependent_properties(
        batch_array, scaled_array, miller.array(ms, data=scales)
    )

    return properties + (batch_data,)


def rmerge_vs_batch(intensities, batches):
    """Determine batches and Rmerge values per batch."""
    assert intensities.size() == batches.size()
    intensities = intensities.map_to_asu()

    merging = intensities.merge_equivalents()
    merged_intensities = merging.array()

    perm = flex.sort_permutation(batches.data())
    batches = batches.data().select(perm)
    intensities = intensities.select(perm)

    pairs = miller.match_multi_indices(
        merged_intensities.indices(), intensities.indices()
    ).pairs()

    def r_merge_per_batch(pairs):
        """Calculate R_merge for the list of (merged-I, I) pairs."""

        merged_indices, unmerged_indices = zip(*pairs)

        unmerged_Ij = intensities.data().select(flex.size_t(unmerged_indices))
        merged_Ij = merged_intensities.data().select(flex.size_t(merged_indices))

        numerator = flex.sum(flex.abs(unmerged_Ij - merged_Ij))
        denominator = flex.sum(unmerged_Ij)

        if denominator > 0:
            return numerator / denominator
        return 0

    return _batch_bins_and_data(batches, pairs, function_to_apply=r_merge_per_batch)


def i_sig_i_vs_batch(intensities, batches):
    """Determine batches and I/sigma values per batch."""
    assert intensities.size() == batches.size()
    assert intensities.sigmas() is not None
    sel = intensities.sigmas() > 0

    i_sig_i = intensities.data().select(sel) / intensities.sigmas().select(sel)
    batches = batches.select(sel)

    perm = flex.sort_permutation(batches.data())
    batches = batches.data().select(perm)
    i_sig_i = i_sig_i.select(perm)

    return _batch_bins_and_data(batches, i_sig_i, function_to_apply=flex.mean)


def scales_vs_batch(scales, batches):
    """Determine batches and scale values per batch."""
    assert scales.size() == batches.size()

    perm = flex.sort_permutation(batches.data())
    batches = batches.data().select(perm)
    scales = scales.data().select(perm)

    return _batch_bins_and_data(batches, scales, function_to_apply=flex.mean)


def _batch_bins_and_data(batches, values, function_to_apply):
    """Apply function to the data from each batch.

    Return the list of the batch bins and the value for each bin.
    """
    batch_bins = []
    data = []

    i_batch_start = 0
    current_batch = flex.min(batches)
    n_ref = batches.size()
    for i_ref in range(n_ref + 1):
        if i_ref == n_ref or batches[i_ref] != current_batch:
            assert batches[i_batch_start:i_ref].all_eq(current_batch)
            values_for_batch = values[i_batch_start:i_ref]
            data.append(function_to_apply(values_for_batch))
            batch_bins.append(current_batch)
            i_batch_start = i_ref
            if i_ref < n_ref:
                current_batch = batches[i_batch_start]
    return batch_bins, data


def make_merging_statistics_summary(dataset_statistics):
    """Format merging statistics information into an output string."""

    # Here use a StringIO to get around excessive padding/whitespace.
    # Avoid using result.show as don't want redundancies printed.
    out = StringIO()

    # First make summary
    make_sub_header("Merging statistics", out=out)
    dataset_statistics.overall.show_summary(out=out)

    # Next make statistics by resolution bin.
    msg = "\n\nStatistics by resolution bin:\n"
    msg += (
        " d_max  d_min   #obs  #uniq   mult.  %comp       <I>  <I/sI>"
        + "    r_mrg   r_meas    r_pim   cc1/2   cc_ano\n"
    )
    for bin_stats in dataset_statistics.bins:
        msg += bin_stats.format() + "\n"
    msg += dataset_statistics.overall.format() + "\n\n"
    out.write(msg)

    # Finally show estimated cutoffs
    dataset_statistics.show_estimated_cutoffs(out=out)
    return out.getvalue()

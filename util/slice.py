from __future__ import absolute_import, division, print_function

import copy

from scitbx.array_family import flex


def slice_experiments(experiments, image_ranges):
    """

    :param experiments
    :type experiments: dxtbx.model.experiment_list.ExperimentList
    :param image_range:
    :type image_range: list of 2-tuples defining scan range for each experiment"""

    # copy the experiments
    experiments = copy.deepcopy(experiments)

    if len(experiments) != len(image_ranges):
        raise Sorry("Input experiment list and image_ranges are not of the same length")

    for exp, sr in zip(experiments, image_ranges):
        if sr is None:
            continue
        im_range = exp.scan.get_image_range()
        if sr[0] < im_range[0] or sr[1] > im_range[1]:
            raise IndexError("requested slice outside current scan range")

        # slicing uses the array range, not the image range
        arr_start = exp.scan.get_array_range()[0]
        beg = sr[0] - 1 - arr_start
        end = sr[1] - arr_start
        exp.scan.swap(exp.scan[beg:end])

    return experiments


def slice_reflections(reflections, image_ranges):
    """

    :param reflections: reflection table of input reflections
    :type reflections: dials.array_family.flex.reflection_table
    :param image_range: list of 2-tuples defining scan range for each experiment
                       id contained within the reflections
    :type image_range: list of 2-tuples defining scan range for each experiment"""

    # copy the reflections
    reflections = copy.deepcopy(reflections)

    to_keep = flex.size_t()
    for iexp, sr in enumerate(image_ranges):

        if sr is None:
            continue
        isel = (reflections["id"] == iexp).iselection()
        frames = (reflections["xyzobs.px.value"].parts()[2]).select(isel)
        # reflns on image n have frames in range [n-1, n)
        in_low_lim = frames >= sr[0] - 1
        in_high_lim = frames < sr[1]
        in_lim = in_low_lim & in_high_lim

        # which indices to keep?
        sub_isel = isel.select(in_lim)
        to_keep.extend(sub_isel)

    # implictly also removes any reflections with ID outside the range of the
    # length of image_ranges
    return reflections.select(to_keep)

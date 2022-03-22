"""
Functions for handling exclusion of images ranges from scans based on a
command line option, as well as obtaining a selection to use for selecting
the corresponding reflections.
"""

from __future__ import annotations

from orderedset import OrderedSet

import iotbx.phil

from dials.array_family import flex

phil_scope = iotbx.phil.parse(
    """
  exclude_images = None
    .type = strings
    .multiple = True
    .help = "Input in the format exp:start:end"
            "Exclude a range of images (start, stop) from the dataset with"
            "experiment identifier exp  (inclusive of frames start, stop)."
            "Multiple ranges can be given in one go, e.g."
            "exclude_images=0:150:200,1:200:250"
            "exclude_images='0:150:200 1:200:250'"
    .short_caption = "Exclude images"
    .expert_level = 1
"""
)


def exclude_image_ranges_from_scans(reflections, experiments, exclude_images):
    """Using an exclude_images phil option, modify the valid image ranges for each
    experiment in the scan (if it exists). Requires experiment identifiers to be
    set already."""
    experiments = set_initial_valid_image_ranges(experiments)
    ranges_to_remove = _parse_exclude_images_commands(
        exclude_images, experiments, reflections
    )
    experiments = _remove_ranges_from_valid_image_ranges(experiments, ranges_to_remove)
    return experiments


def get_valid_image_ranges(experiments):
    """Extract valid image ranges from experiments, returning None if no scan"""
    valid_images_ranges = []
    for exp in experiments:
        if exp.scan:
            valid_images_ranges.append(exp.scan.get_valid_image_ranges(exp.identifier))
        else:
            valid_images_ranges.append(None)
    return valid_images_ranges


def set_initial_valid_image_ranges(experiments):
    """Set the valid_image_range for each experiment to be the scan image range.
    Kept separate from dxtbx.scan object as requires experiment indentifiers.
    Also this function can be called for a mix of sequences and scanless experiments.
    """
    for exp in experiments:
        if exp.scan:
            if not exp.scan.get_valid_image_ranges(exp.identifier):
                exp.scan.set_valid_image_ranges(
                    exp.identifier, [exp.scan.get_image_range()]
                )
    return experiments


def get_selection_for_valid_image_ranges(reflection_table, experiment):
    """Determine a selection for the reflection table corresponding to reflections
    that are located in valid image ranges (according to zobs.px.value)."""
    if experiment.scan:
        valid_ranges = experiment.scan.get_valid_image_ranges(experiment.identifier)
        if valid_ranges:
            valid_mask = flex.bool(reflection_table.size(), False)
            z = reflection_table["xyzobs.px.value"].parts()[2]
            for valid_range in valid_ranges:  # i in range(0, len(valid_ranges), 2):
                mask1 = z >= valid_range[0] - 1.0  # -1.0 as image valid_ranges[i] \
                # has z is valid_ranges[i] - 1 to valid_ranges[i]
                mask2 = z <= valid_range[1]
                valid_mask.set_selected(mask1 & mask2, True)
            return valid_mask
    return flex.bool(reflection_table.size(), True)  # else say all valid


def _parse_exclude_images_commands(commands, experiments, reflections):
    """Parse a list of list of command line options.

    e.g. commands = [['1:101:200'], ['0:201:300']]
    or commands = [[101:200]] allowable for a single experiment.
    builds and returns a list of tuples (exp_id, (start, stop))


    """
    ranges_to_remove = []
    for com in commands:
        if len(com) > 1:
            sub = com
        elif "," in com[0]:
            sub = com[0].split(",")
        else:
            sub = com
        for subvals in sub:
            vals = subvals.split(":")
            if len(vals) == 2:
                if len(experiments) > 1:
                    raise ValueError(
                        "Exclude images must be in the form experimentnumber:start:stop for multiple experiments"
                    )
                else:
                    ranges_to_remove.append(
                        (experiments[0].identifier, (int(vals[0]), int(vals[1])))
                    )
            else:
                if len(vals) != 3:
                    raise ValueError(
                        "Exclude images must be input in the form experimentnumber:start:stop, or start:stop for a single experiment"
                        + "Multiple ranges can be specified by comma or space separated values e.g 0:100:150,1:120:200"
                    )
                dataset_id = int(vals[0])
                for table in reflections:
                    if dataset_id in table.experiment_identifiers():
                        expid = table.experiment_identifiers()[dataset_id]
                        ranges_to_remove.append((expid, (int(vals[1]), int(vals[2]))))
                        break
    return ranges_to_remove


def _remove_ranges_from_valid_image_ranges(experiments, ranges_to_remove):
    """Update the valid_image_ranges for each experiment in the scans. Uses set
    arithmetic to determine image ranges to keep and sets a new list of ranges in
    the scan.valid_image_ranges dict."""
    ids = list(experiments.identifiers())
    for r in ranges_to_remove:
        idx = ids.index(r[0])
        exp = experiments[idx]
        if not exp.scan:
            raise ValueError("Trying to exclude a scanless experiment")
        current_range = exp.scan.get_valid_image_ranges(
            exp.identifier
        )  # list of tuples
        # use set arithmetic on image numbers to work out images to keep
        exclude = OrderedSet(range(r[1][0], r[1][1] + 1))
        current_sets = [
            OrderedSet(range(current_range[i][0], current_range[i][1] + 1))
            for i in range(0, len(current_range))
        ]
        current_images = OrderedSet()
        for s in current_sets:
            current_images = current_images | s
        new_valid_images = list(current_images - exclude)
        # now convert images to ranges for storing in the scan
        if new_valid_images:
            new_valid_ranges = [new_valid_images[0]]
            previous = new_valid_images[0]
            for i in new_valid_images[1:]:
                if i - previous > 1:
                    new_valid_ranges.append(previous)
                    new_valid_ranges.append(i)
                previous = i
            # add end one
            new_valid_ranges.append(new_valid_images[-1])
            assert len(new_valid_ranges) % 2 == 0
            # convert to list of tuples
            valid_ranges = [
                (new_valid_ranges[i], new_valid_ranges[i + 1])
                for i in range(0, len(new_valid_ranges), 2)
            ]
        else:
            valid_ranges = []
        exp.scan.set_valid_image_ranges(exp.identifier, valid_ranges)
    return experiments


"""
Functions for scaling
"""


def exclude_image_ranges_for_scaling(reflections, experiments, exclude_images):
    """Set the initial valid ranges, then exclude in the exp.scan and set
    user_excluded_in_scaling flags."""
    experiments = exclude_image_ranges_from_scans(
        reflections, experiments, exclude_images
    )
    for refl, exp in zip(reflections, experiments):
        sel = get_selection_for_valid_image_ranges(refl, exp)
        refl.set_flags(~sel, refl.flags.user_excluded_in_scaling)
        refl.unset_flags(~sel, refl.flags.scaled)
    return reflections, experiments

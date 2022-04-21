"""
Functions to help with calculating batch properties for experiments objects.
"""

from __future__ import annotations

import logging

from dials.array_family import flex

logger = logging.getLogger("dials")


class batch_manager:
    def __init__(self, batches, batch_params):
        # batch params is a list of dicts with "id" and "range" - used to be
        # a 'scope extract' object
        self.batch_params = sorted(batch_params, key=lambda b: b["range"][0])
        self.batches = batches
        self.reduced_batches, self._batch_increments = self._reduce()

    def _reduce(self):
        reduced_batches = flex.int(self.batches)
        batch_increments = []
        incr = 0
        for batch in self.batch_params:
            sel = (reduced_batches >= batch["range"][0]) & (
                reduced_batches <= batch["range"][1]
            )
            reduced_batches.set_selected(
                sel, reduced_batches.select(sel) - (batch["range"][0] - incr) + 1
            )
            batch_increments.append(incr)
            incr += batch["range"][1] - batch["range"][0] + 1
        assert len(set(reduced_batches)) == len(reduced_batches)
        return list(reduced_batches), batch_increments

    def batch_plot_shapes_and_annotations(self):
        light_grey = "#d3d3d3"
        grey = "#808080"
        shapes = []
        annotations = []
        batches = flex.int(self.batches)
        text = flex.std_string(batches.size())
        for i, batch in enumerate(self.batch_params):
            fillcolor = [light_grey, grey][i % 2]  # alternate colours
            shapes.append(
                {
                    "type": "rect",
                    # x-reference is assigned to the x-values
                    "xref": "x",
                    # y-reference is assigned to the plot paper [0,1]
                    "yref": "paper",
                    "x0": self._batch_increments[i],
                    "y0": 0,
                    "x1": self._batch_increments[i]
                    + (batch["range"][1] - batch["range"][0]),
                    "y1": 1,
                    "fillcolor": fillcolor,
                    "opacity": 0.2,
                    "line": {"width": 0},
                }
            )
            annotations.append(
                {
                    # x-reference is assigned to the x-values
                    "xref": "x",
                    # y-reference is assigned to the plot paper [0,1]
                    "yref": "paper",
                    "x": self._batch_increments[i]
                    + (batch["range"][1] - batch["range"][0]) / 2,
                    "y": 1,
                    "text": f"{batch['id']}",
                    "showarrow": False,
                    "yshift": 20,
                    # 'arrowhead': 7,
                    # 'ax': 0,
                    # 'ay': -40
                }
            )
            sel = (batches >= batch["range"][0]) & (batches <= batch["range"][1])
            text.set_selected(
                sel,
                flex.std_string(
                    [
                        f"{batch['id']}: {j - batch['range'][0] + 1}"
                        for j in batches.select(sel)
                    ]
                ),
            )
        return shapes, annotations, list(text)


def assign_batches_to_reflections(reflections, batch_offsets):
    """Assign a 'batch' column to the reflection table"""
    for batch_offset, refl in zip(batch_offsets, reflections):
        xdet, ydet, zdet = [flex.double(x) for x in refl["xyzobs.px.value"].parts()]
        # compute BATCH values - floor() to get (fortran) image captured within
        #                        +1     because FORTRAN counting; zdet+1=image_index
        #                        +off   because            image_index+o=batch
        refl["batch"] = (flex.floor(zdet).iround() + 1) + batch_offset
    return reflections


def get_batch_ranges(experiments, batch_offsets):
    """Get batch ranges for a list of experiments and offsets"""
    batch_ranges = []
    assert len(experiments) == len(batch_offsets)
    image_ranges = get_image_ranges(experiments)
    for batch_offset, image_range in zip(batch_offsets, image_ranges):
        batch_ranges.append(
            (image_range[0] + batch_offset, image_range[1] + batch_offset)
        )
    return batch_ranges


def get_image_ranges(experiments):
    """Get image ranges for a list of experiments (including scanless exp.)"""
    # Note, if set to 1,1,for scanless experiments then first batch offset in
    # _calculate_batch_offsets is zero below, bad!
    return [e.scan.get_image_range() if e.scan else (0, 0) for e in experiments]


def calculate_batch_offsets(experiment_list):
    """Take a list of experiments and resolve and return the batch offsets.
    First adds an image_range property as not all experiments have scans."""
    image_ranges = get_image_ranges(experiment_list)
    offsets = _calculate_batch_offsets(image_ranges)
    return offsets


def set_batch_offsets(experiment_list, batch_offsets):
    """Set batch offsets in scan objects. Don't need to set anything for
    scanless experiments, as these are not used with the batch system."""
    for exp, offset in zip(experiment_list, batch_offsets):
        if exp.scan:
            exp.scan.set_batch_offset(offset)


def _calculate_batch_offsets(image_ranges):
    """Take a list of (modified) experiments and resolve and return the batch
    offsets.

    This is the number added to the image number to give the
    batch number, such that:
    - Each experiment has a unique, nonoverlapping, nonconsecutive range
    - None are zero
    - Image number ranges are kept if at all possible
    """

    experiments_to_shift = []
    existing_ranges = set()
    maximum_batch_number = 0
    batch_offsets = [0] * len(image_ranges)

    # Handle zeroth shifts and kept ranges
    for i, image_range in enumerate(image_ranges):
        ilow, ihigh = image_range
        # Check assumptions
        assert ilow <= ihigh, "Inverted image order!?"
        assert ilow >= 0, "Negative image indices are not expected"
        # Don't emit zero: Causes problems with C/fortran number conversion
        if ilow == 0:
            ilow, ihigh = ilow + 1, ihigh + 1
        # If we overlap with anything, then process later
        if any(ilow < high + 1 and ihigh >= low - 1 for low, high in existing_ranges):
            experiments_to_shift.append((i, image_range))
        else:
            batch_offsets[i] = ilow - image_range[0]
            existing_ranges.add((ilow, ihigh))
            maximum_batch_number = max(maximum_batch_number, ihigh)
    # Now handle all the experiments that overlapped by pushing them higher
    for i, image_range in experiments_to_shift:
        start_number = _next_epoch(maximum_batch_number)
        range_width = image_range[1] - image_range[0] + 1
        end_number = start_number + range_width - 1
        batch_offsets[i] = start_number - image_range[0]
        maximum_batch_number = end_number
    return batch_offsets


def _next_epoch(val):
    """Find the next number above the existing value that ends in 1, that is
    not consecutive with the current value."""
    if val % 100 == 99:
        return val + 2
    elif val % 100 == 0:
        return val + 101
    else:
        rem = val % 100
        return val - rem + 101

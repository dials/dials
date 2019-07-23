#!/usr/bin/env python
#
# dials.slice_sweep.py
#
#  Copyright (C) 2014 Diamond Light Source and STFC Rutherford Appleton
#  Laboratory, UK.
#
#  Author: David Waterman
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import absolute_import, division, print_function

from os.path import basename, splitext

from dials.algorithms.refinement.refinement_helpers import calculate_frame_numbers
from dxtbx.model.experiment_list import ExperimentList
from dials.util import Sorry
from dials.util.slice import slice_experiments, slice_reflections

help_message = """

Slice a sweep to produce a smaller sweep within the bounds of the original. If
experiments or experiments are provided, modify the scan objects within these. If
reflections are provided, remove reflections outside the provided image ranges.
Each image_range parameter refers to a single experiment ID, counting up from
zero. Any reflections with experiment ID not matched by a image_range parameter
are removed.

Examples::

  dials.slice_sweep models.expt observations.refl "image_range=1 20"

  dials.slice_sweep models.expt "image_range=1 20"

  # two experiments and reflections with IDs '0' and '1'
  dials.slice_sweep models.expt observations.refl \
    "image_range=1 20" "image_range=5 30"

"""

# The phil scope
from libtbx.phil import parse

phil_scope = parse(
    """

  output {

    reflections_filename = None
      .type = str
      .help = "The filename for output reflections sliced to remove those"
              "outside the reduced image range. By default generated"
              "automatically from the input name"

    experiments_filename = None
      .type = str
      .help = "The filename for the output experiments with sliced scans.
               By default generated automatically from the input name"

  }

  image_range = None
    .help = "Range in images to slice a sweep. The number of arguments"
            "must be a factor of two. Each pair of arguments gives a range"
            "that follows C conventions (e.g. j0 <= j < j1) when slicing the"
            "reflections by observed centroid."
    .type = ints(size=2)
    .multiple = True

  block_size = None
    .type = float
    .help = "Overrides image_range if present. This option splits each sweep"
            "into the nearest integer number of equal size blocks close to"
            "block_size degrees in width"

"""
)


def calculate_block_ranges(scan, block_size):
    """

    :param scans
    :type a scan object
    :param block_size:
    :type block_size: target block size in degrees"""

    image_ranges = []
    nimages = scan.get_num_images()
    osc_range = scan.get_oscillation_range(deg=True)
    osc_width = abs(osc_range[1] - osc_range[0])
    nblocks = max(int(round(osc_width / block_size)), 1)
    nblocks = min(nblocks, nimages)
    # equal sized blocks except the last one that may contain extra images
    # to make up the remainder
    nimages_per_block = [nimages // nblocks] * (nblocks - 1) + [
        nimages // nblocks + nimages % nblocks
    ]
    start = scan.get_image_range()[0]
    for nim in nimages_per_block:
        image_ranges.append((start, start + nim - 1))
        start += nim

    return image_ranges


class Script(object):
    """A class for running the script."""

    def __init__(self):
        """Initialise the script."""
        from dials.util.options import OptionParser

        usage = (
            "usage: dials.slice_sweep [options] [param.phil] "
            "models.expt observations.refl"
        )

        # Create the parser
        self.parser = OptionParser(
            usage=usage,
            phil=phil_scope,
            read_reflections=True,
            read_experiments=True,
            check_format=False,
            epilog=help_message,
        )

    def run(self):
        """Execute the script."""

        from dials.util.options import flatten_reflections, flatten_experiments

        # Parse the command line
        params, options = self.parser.parse_args(show_diff_phil=True)
        reflections = flatten_reflections(params.input.reflections)
        experiments = flatten_experiments(params.input.experiments)

        # Try to load the models and data
        slice_exps = len(experiments) > 0
        slice_refs = len(reflections) > 0

        # Catch case of nothing to do
        if not any([slice_exps, slice_refs]):
            print("No suitable input provided")
            self.parser.print_help()
            return

        if reflections:
            if len(reflections) > 1:
                raise Sorry("Only one reflections list can be imported at present")
            reflections = reflections[0]

            # calculate frame numbers if needed
            if experiments:
                reflections = calculate_frame_numbers(reflections, experiments)

            # if we still don't have the right column give up
            if "xyzobs.px.value" not in reflections:
                raise Sorry(
                    "These reflections do not have frame numbers set, and "
                    "there are no experiments provided to calculate these."
                )

        # set trivial case where no scan range is provided at all
        if not params.image_range:
            params.image_range = [None]

        # check if slicing into blocks
        if params.block_size is not None:
            if slice_exps:
                if len(experiments) > 1:
                    raise Sorry(
                        "For slicing into blocks please provide a single " "scan only"
                    )
                scan = experiments[0].scan

            # Having extracted the scan, calculate the blocks
            params.image_range = calculate_block_ranges(scan, params.block_size)

            # Do the slicing then recombine
            if slice_exps:
                sliced = [
                    slice_experiments(experiments, [sr])[0] for sr in params.image_range
                ]
                sliced_experiments = ExperimentList()
                for exp in sliced:
                    sliced_experiments.append(exp)

            # slice reflections if present
            if slice_refs:
                sliced = [
                    slice_reflections(reflections, [sr]) for sr in params.image_range
                ]
                sliced_reflections = sliced[0]
                for i, rt in enumerate(sliced[1:]):
                    rt["id"] += i + 1  # set id
                    sliced_reflections.extend(rt)

        else:
            # slice each dataset into the requested subset
            if slice_exps:
                sliced_experiments = slice_experiments(experiments, params.image_range)
            if slice_refs:
                sliced_reflections = slice_reflections(reflections, params.image_range)

        # Save sliced experiments
        if slice_exps:
            output_experiments_filename = params.output.experiments_filename
            if output_experiments_filename is None:
                # take first filename as template
                bname = basename(params.input.experiments[0].filename)
                bname = splitext(bname)[0]
                if not bname:
                    bname = "experiments"
                if len(params.image_range) == 1 and params.image_range[0] is not None:
                    ext = "_{0}_{1}.expt".format(*params.image_range[0])
                else:
                    ext = "_sliced.expt"
                output_experiments_filename = bname + ext
            print("Saving sliced experiments to {}".format(output_experiments_filename))

            from dxtbx.model.experiment_list import ExperimentListDumper

            dump = ExperimentListDumper(sliced_experiments)
            dump.as_json(output_experiments_filename)

        # Save sliced reflections
        if slice_refs:
            output_reflections_filename = params.output.reflections_filename
            if output_reflections_filename is None:
                # take first filename as template
                bname = basename(params.input.reflections[0].filename)
                bname = splitext(bname)[0]
                if not bname:
                    bname = "reflections"
                if len(params.image_range) == 1 and params.image_range[0] is not None:
                    ext = "_{0}_{1}.refl".format(*params.image_range[0])
                else:
                    ext = "_sliced.refl"
                output_reflections_filename = bname + ext

            print(
                "Saving sliced reflections to {0}".format(output_reflections_filename)
            )
            sliced_reflections.as_pickle(output_reflections_filename)

        return


if __name__ == "__main__":
    from dials.util import halraiser

    try:
        script = Script()
        script.run()
    except Exception as e:
        halraiser(e)

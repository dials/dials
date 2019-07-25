#!/usr/bin/env python
#
# dials.compute_delta_cchalf
#
#  Copyright (C) 2018 Diamond Light Source
#
#  Authors: James Parkhurst, James Beilsten-Edmands
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#
from __future__ import absolute_import, division, print_function
import collections
import logging
from math import sqrt
import matplotlib
from iotbx.reflection_file_reader import any_reflection_file
from libtbx.phil import parse
import dials.util
from dials.algorithms.statistics.delta_cchalf import PerImageCChalfStatistics
from dials.array_family import flex
from dials.util import Sorry
from dials.util.exclude_images import exclude_image_ranges_for_scaling
from dials.util.multi_dataset_handling import select_datasets_on_ids

matplotlib.use("Agg")
from matplotlib import pylab

logger = logging.getLogger("dials.command_line.compute_delta_cchalf")

help_message = """

This program computes the delta cchalf excluding images

"""

# Set the phil scope
phil_scope = parse(
    """

  input {

    mtzfile = None
      .type = str
      .help = "We can also import an MTZ file"

  }

  mode = *dataset image_group
    .type = choice
    .help = "Perform analysis on whole datasets or batch groups"

  group_size = 10
    .type = int(value_min=1)
    .help = "The number of images to group together when calculating delta"
            "cchalf in image_group mode"

  output {

    experiments = "filtered.expt"
      .type = str
      .help = "The filtered experiments file"

    reflections = "filtered.refl"
      .type = str
      .help = "The filtered reflections file"

    table = "delta_cchalf.dat"
      .type = str
      .help = "A file with delta cchalf values"
  }

  nbins = 10
    .type = int(value_min=1)
    .help = "The number of resolution bins to use"

  dmin = None
    .type = float
    .help = "The maximum resolution"

  dmax = None
    .type = float
    .help = "The minimum resolution"

  stdcutoff = 4.0
    .type = float
    .help = "Datasets with a delta cc half below (mean - stdcutoff*std) are removed"

  output {

    log = 'dials.compute_delta_cchalf.log'
      .type = str
      .help = "The log filename"

    debug_log = 'dials.compute_delta_cchalf.debug.log'
      .type = str
      .help = "The debug log filename"
  }
"""
)


class Script(object):
    """A class for running the script."""

    def __init__(self, params, experiments, reflections):
        """Initialise the script."""
        self.experiments = experiments
        self.reflections = reflections
        self.params = params
        self.delta_cchalf_i = {}
        self.results_summary = {
            "dataset_removal": {
                "mode": self.params.mode,
                "stdcutoff": self.params.stdcutoff,
            }
        }
        # Set up a named tuple
        self.DataRecord = collections.namedtuple(
            "DataRecord",
            (
                "unit_cell",
                "space_group",
                "miller_index",
                "dataset",
                "intensity",
                "variance",
                "identifiers",
                "images",
            ),
        )

    def prepare_data(self):
        """Prepare the data into a DataRecord."""
        if self.params.mode == "image_group":
            for exp in self.experiments:
                if not exp.scan:
                    raise Sorry("Cannot use mode=image_group with scanless experiments")
        if len(self.experiments) > 0:
            if len(self.experiments) > 0 and len(self.reflections) == 1:
                data = self.read_experiments(self.experiments, self.reflections[0])
            elif len(self.experiments) > 0 and len(self.experiments) == len(
                self.reflections
            ):
                # need to join together reflections
                joint_table = flex.reflection_table()
                for table in self.reflections:
                    joint_table.extend(table)
                self.reflections = [joint_table]
                data = self.read_experiments(self.experiments, self.reflections[0])
            else:
                raise Sorry("Unequal number of reflection tables and experiments")
        elif self.params.input.mtzfile is not None:
            data = self.read_mtzfile(self.params.input.mtzfile)
        else:
            return SystemExit
        return data

    def run(self):
        """Run the delta_cc_half algorithm."""

        data = self.prepare_data()

        # Create the statistics object
        statistics = PerImageCChalfStatistics(
            data.miller_index,
            data.identifiers,
            data.dataset,
            data.images,
            data.intensity,
            data.variance,
            data.unit_cell,
            data.space_group,
            self.params.nbins,
            self.params.dmin,
            self.params.dmax,
            self.params.mode,
            self.params.group_size,
        )

        self.delta_cchalf_i = statistics.delta_cchalf_i()
        self.results_summary["mean_cc_half"] = statistics._cchalf_mean
        # Print out the datasets in order of delta cc 1/2
        sorted_datasets, sorted_cc_half_values = self.sort_deltacchalf_values(
            self.delta_cchalf_i, self.results_summary
        )

        # Write a text file with delta cchalf value
        logger.info("Writing table to %s", self.params.output.table)
        with open(self.params.output.table, "w") as outfile:
            for dataset, cchalf in zip(sorted_datasets, sorted_cc_half_values):
                outfile.write("%d %f\n" % (dataset, cchalf))

        # Remove datasets based on delta cc1/2
        if self.experiments and len(self.reflections) == 1:
            cutoff_value = self._calculate_cutoff_value(
                self.delta_cchalf_i, self.params.stdcutoff
            )
            self.results_summary["dataset_removal"].update(
                {"cutoff_value": cutoff_value}
            )

            below_cutoff = sorted_cc_half_values < cutoff_value
            ids_to_remove = sorted_datasets.select(below_cutoff)

            if self.params.mode == "dataset":
                filtered_reflections = self.remove_datasets_below_cutoff(
                    self.experiments,
                    self.reflections[0],
                    ids_to_remove,
                    self.results_summary,
                )
            elif self.params.mode == "image_group":
                filtered_reflections = self.remove_image_ranges_below_cutoff(
                    self.experiments,
                    self.reflections[0],
                    ids_to_remove,
                    statistics.image_group_to_expid_and_range,
                    statistics.expid_to_image_groups,
                    self.results_summary,
                )
            self.reflections = [filtered_reflections]

    def read_experiments(self, experiments, reflections):
        """
        Get information from experiments and reflections

        """

        # Get space group and unit cell
        space_group = None
        unit_cell = []
        exp_identifiers = []
        for e in experiments:
            if space_group is None:
                space_group = e.crystal.get_space_group()
            else:
                assert (
                    space_group.type().number()
                    == e.crystal.get_space_group().type().number()
                )
            unit_cell.append(e.crystal.get_unit_cell())
            exp_identifiers.append(e.identifier)
        # get a list of the ids from the reflection table corresponding to exp_ids
        identifiers = []
        for expit in exp_identifiers:
            for k in reflections.experiment_identifiers().keys():
                if reflections.experiment_identifiers()[k] == expit:
                    identifiers.append(k)
                    break

        # Selection of reflections
        selection = ~(
            reflections.get_flags(reflections.flags.bad_for_scaling, all=False)
        )
        outliers = reflections.get_flags(reflections.flags.outlier_in_scaling)
        reflections = reflections.select(selection & ~outliers)

        # Scale factor
        inv_scale_factor = reflections["inverse_scale_factor"]
        selection = inv_scale_factor > 0
        reflections = reflections.select(selection)
        inv_scale_factor = reflections["inverse_scale_factor"]

        # Get the reflection data
        index = reflections["id"]

        miller_index = reflections["miller_index"]
        intensity = reflections["intensity.scale.value"] / inv_scale_factor
        variance = reflections["intensity.scale.variance"] / inv_scale_factor ** 2
        # calculate image number of observation (e.g 0.0 <= z < 1.0), image = 1
        images = flex.floor(reflections["xyzobs.px.value"].parts()[2]).iround() + 1
        # Get the MTZ file
        return self.DataRecord(
            unit_cell=unit_cell,
            space_group=space_group,
            miller_index=miller_index,
            dataset=index,
            intensity=intensity,
            variance=variance,
            identifiers=identifiers,
            images=images,
        )

    def read_mtzfile(self, filename):
        """
        Read the mtz file

        """
        # Read the mtz file
        reader = any_reflection_file(filename)

        # Get the columns as miller arrays
        miller_arrays = reader.as_miller_arrays(merge_equivalents=False)

        # Select the desired columns
        intensities = None
        batches = None
        for array in miller_arrays:
            if array.info().labels == ["I", "SIGI"]:
                intensities = array
            if array.info().labels == ["BATCH"]:
                batches = array
        assert intensities is not None
        assert batches is not None
        assert len(batches.data()) == len(intensities.data())

        # Get the unit cell and space group
        unit_cell = intensities.unit_cell()
        space_group = intensities.crystal_symmetry().space_group()

        # The reflection data
        miller_index = intensities.indices()
        batch = batches.data()
        intensity = intensities.data()
        variance = intensities.sigmas() ** 2

        # Create unit cell list
        min_batch = min(batch)
        dataset = batch - min_batch
        num_datasets = max(dataset) + 1
        unit_cell_list = [unit_cell for _ in range(num_datasets)]

        # Get the MTZ file
        return self.DataRecord(
            unit_cell=unit_cell_list,
            space_group=space_group,
            miller_index=miller_index,
            dataset=dataset,
            intensity=intensity,
            variance=variance,
        )

    @staticmethod
    def sort_deltacchalf_values(delta_cchalf_i, results_summary):
        """Return the sorted datasets and cchalf values.

        Also add the sorted lists to the results summary. Datasets are sorted
        from low to high based on deltacchalf values."""
        datasets = list(delta_cchalf_i.keys())
        sorted_index = sorted(
            range(len(datasets)), key=lambda x: delta_cchalf_i[datasets[x]]
        )

        # sorted by deltacchalf from low to high
        sorted_cc_half_values = flex.double([])
        sorted_datasets = flex.int([])
        for i in sorted_index:
            val = delta_cchalf_i[datasets[i]]
            logger.info("Dataset: %d, Delta CC 1/2: %.3f", datasets[i], 100 * val)
            sorted_cc_half_values.append(val)
            sorted_datasets.append(datasets[i])

        results_summary["per_dataset_delta_cc_half_values"] = {
            "datasets": list(sorted_datasets),
            "delta_cc_half_values": list(sorted_cc_half_values),
        }

        return sorted_datasets, sorted_cc_half_values

    @staticmethod
    def _calculate_cutoff_value(delta_cchalf_i, stdcutoff):
        Y = list(delta_cchalf_i.values())
        mean = sum(Y) / len(Y)
        sdev = sqrt(sum((yy - mean) ** 2 for yy in Y) / len(Y))
        logger.info("\nmean delta_cc_half %s", (mean * 100))
        logger.info("stddev delta_cc_half %s", (sdev * 100))
        cutoff_value = mean - stdcutoff * sdev
        logger.info("cutoff value: %s \n", (cutoff_value * 100))
        return cutoff_value

    @staticmethod
    def remove_image_ranges_below_cutoff(
        experiments,
        reflections,
        ids_to_remove,
        image_group_to_expid_and_range,
        expid_to_image_groups,
        results_summary,
    ):
        """Remove image ranges from the datasets."""
        n_valid_reflections = reflections.get_flags(
            reflections.flags.bad_for_scaling, all=False
        ).count(False)

        experiments_to_delete = []
        exclude_images = []
        image_ranges_removed = []  # track for results summary
        n_removed_this_cycle = 1
        while n_removed_this_cycle != 0:
            other_potential_ids_to_remove = []
            n_removed_this_cycle = 0
            for id_ in sorted(ids_to_remove):
                exp_id, image_range = image_group_to_expid_and_range[
                    id_
                ]  # numerical id
                identifier = reflections.experiment_identifiers()[exp_id]
                if expid_to_image_groups[exp_id][-1] == id_:  # is last group
                    image_ranges_removed.append([image_range, exp_id])
                    logger.info(
                        "Removing image range %s from experiment %s",
                        image_range,
                        identifier,
                    )
                    exclude_images.append(
                        [
                            identifier
                            + ":"
                            + str(image_range[0])
                            + ":"
                            + str(image_range[1])
                        ]
                    )
                    del expid_to_image_groups[exp_id][-1]
                    n_removed_this_cycle += 1
                else:
                    other_potential_ids_to_remove.append(id_)
            ids_to_remove = other_potential_ids_to_remove
        for id_ in other_potential_ids_to_remove:
            exp_id, image_range = image_group_to_expid_and_range[id_]
            identifier = reflections.experiment_identifiers()[exp_id]
            logger.info(
                """Image range %s from experiment %s is below the cutoff, but not at the end of a sweep.""",
                image_range,
                identifier,
            )

        # Now remove individual batches
        if -1 in reflections["id"]:
            reflections = reflections.select(reflections["id"] != -1)
        reflection_list = reflections.split_by_experiment_id()
        reflection_list, experiments = exclude_image_ranges_for_scaling(
            reflection_list, experiments, exclude_images
        )
        # if a whole experiment has been excluded: need to remove it here

        for exp in experiments:
            if not exp.scan.get_valid_image_ranges(
                exp.identifier
            ):  # if all removed above
                experiments_to_delete.append(exp.identifier)
        if experiments_to_delete:
            experiments, reflection_list = select_datasets_on_ids(
                experiments, reflection_list, exclude_datasets=experiments_to_delete
            )
        assert len(reflection_list) == len(experiments)

        output_reflections = flex.reflection_table()
        for r in reflection_list:
            output_reflections.extend(r)

        n_valid_filtered_reflections = output_reflections.get_flags(
            output_reflections.flags.bad_for_scaling, all=False
        ).count(False)
        results_summary["dataset_removal"].update(
            {
                "image_ranges_removed": image_ranges_removed,
                "experiments_fully_removed": experiments_to_delete,
                "n_reflections_removed": n_valid_reflections
                - n_valid_filtered_reflections,
            }
        )
        return output_reflections

    @staticmethod
    def remove_datasets_below_cutoff(
        experiments, reflections, ids_to_remove, results_summary
    ):
        """Remove the datasets with ids in ids_to_remove.

        Remove from the experiemnts and reflections and add information to the
        results summary dict.

        Returns:
          output_reflections: The reflection table with data removed.
        """
        n_valid_reflections = reflections.get_flags(
            reflections.flags.bad_for_scaling, all=False
        ).count(False)

        datasets_to_remove = []
        for id_ in sorted(ids_to_remove):
            logger.info("Removing dataset %d", id_)
            datasets_to_remove.append(reflections.experiment_identifiers()[id_])
        output_reflections = reflections.remove_on_experiment_identifiers(
            datasets_to_remove
        )
        experiments.remove_on_experiment_identifiers(datasets_to_remove)
        output_reflections.assert_experiment_identifiers_are_consistent(experiments)

        n_valid_filtered_reflections = output_reflections.get_flags(
            output_reflections.flags.bad_for_scaling, all=False
        ).count(False)
        results_summary["dataset_removal"].update(
            {
                "experiments_fully_removed": datasets_to_remove,
                "n_reflections_removed": n_valid_reflections
                - n_valid_filtered_reflections,
            }
        )
        return output_reflections

    def write_experiments_and_reflections(self):
        """Save the reflections and experiments data."""
        if self.experiments and len(self.reflections) == 1:
            logger.info(
                "Saving %d reflections to %s",
                len(self.reflections[0]),
                self.params.output.reflections,
            )
            self.reflections[0].as_pickle(self.params.output.reflections)
            from dxtbx.model.experiment_list import ExperimentListDumper

            logger.info("Saving the experiments to %s", self.params.output.experiments)
            dump = ExperimentListDumper(self.experiments)
            with open(self.params.output.experiments, "w") as outfile:
                outfile.write(dump.as_json())

    def plot_data(self):
        """Plot histogram and line plot of cc half values."""
        fig, ax = pylab.subplots()
        ax.hist(list(self.delta_cchalf_i.values()))
        ax.set_xlabel("Delta CC 1/2")
        fig.savefig("plot1.png")

        X = list(self.delta_cchalf_i.keys())
        Y = list(self.delta_cchalf_i.values())
        fig, ax = pylab.subplots()
        ax.plot(X, Y)
        ax.set_xlabel("Dataset number")
        ax.set_ylabel("Delta CC 1/2")
        fig.savefig("plot2.png")


def run(args=None, phil=phil_scope):
    """Run the command-line script."""
    import dials.util.log
    from dials.util.options import OptionParser
    from dials.util.options import flatten_reflections
    from dials.util.options import flatten_experiments

    usage = "dials.compute_delta_cchalf [options] scaled.expt scaled.refl"

    parser = OptionParser(
        usage=usage,
        phil=phil,
        epilog=__doc__,
        read_experiments=True,
        read_reflections=True,
        check_format=False,
    )

    params, _ = parser.parse_args(args=args, show_diff_phil=False)

    dials.util.log.config(info=params.output.log, debug=params.output.debug_log)

    experiments = flatten_experiments(params.input.experiments)
    reflections = flatten_reflections(params.input.reflections)

    script = Script(params, experiments, reflections)
    script.run()
    script.write_experiments_and_reflections()
    script.plot_data()


if __name__ == "__main__":
    with dials.util.show_mail_on_error():
        run()

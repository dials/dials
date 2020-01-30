from __future__ import absolute_import, division, print_function

import logging
from collections import defaultdict
from math import sqrt
import matplotlib
from iotbx import mtz
from libtbx.phil import parse
import dials.util
from cctbx import uctbx
from dials.algorithms.statistics.delta_cchalf import PerImageCChalfStatistics
from dials.array_family import flex
from dials.util import Sorry
from dials.util.exclude_images import exclude_image_ranges_for_scaling
from dials.util.multi_dataset_handling import select_datasets_on_identifiers
from dials.util.filter_reflections import filter_reflection_table

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

  mtz {
    batch_offset = None
      .type = int
      .help = "The batch offset between consecutive datasets in the mtz file."
  }

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
        self.group_to_datasetid_and_range = {}
        self.datasetid_to_groups = defaultdict(list)

    def prepare_data(self):
        """Prepare the data into a DataRecord."""
        if self.params.mode == "image_group":
            for exp in self.experiments:
                if not exp.scan:
                    raise Sorry("Cannot use mode=image_group with scanless experiments")
        if self.experiments:
            if self.experiments and len(self.reflections) == 1:
                data = self.read_experiments(self.experiments, self.reflections[0])
            elif self.experiments and len(self.experiments) == len(self.reflections):
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
            raise SystemExit
        return data

    def run(self):
        """Run the delta_cc_half algorithm."""

        table, unit_cell_list, space_group = self.prepare_data()

        if len(unit_cell_list) > 1:
            # calc mean unit cell
            mean_parameters = [0, 0, 0, 0, 0, 0]
            for uc in unit_cell_list:
                for i in range(6):
                    mean_parameters[i] += uc.parameters()[i]
            for i in range(6):
                mean_parameters[i] /= len(unit_cell_list)
            mean_unit_cell = uctbx.unit_cell(mean_parameters)
        else:
            mean_unit_cell = unit_cell_list[0]

        if self.params.mode == "dataset":
            table["group"] = table["dataset"]
        elif self.params.mode == "image_group":
            # set up tracking of groups to expts

            image_groups = flex.int(table["dataset"].size(), 0)
            counter = 0
            for id_ in set(table["dataset"]):
                sel = table["dataset"] == id_
                images_in_dataset = table["image"].select(sel)
                unique_images = set(images_in_dataset)
                min_img, max_img = (min(unique_images), max(unique_images))
                group_starts = list(range(min_img, max_img + 1, self.params.group_size))
                if len(group_starts) > 1:
                    if max_img - group_starts[-1] < (self.params.group_size - 1):
                        del group_starts[-1]

                for i, start in enumerate(group_starts[:-1]):
                    group_sel = (images_in_dataset >= start) & (
                        images_in_dataset < group_starts[i + 1]
                    )
                    image_groups.set_selected(
                        (sel.iselection().select(group_sel)), counter
                    )
                    self.group_to_datasetid_and_range[counter] = (
                        id_,
                        (start, group_starts[i + 1] - 1),
                    )
                    self.datasetid_to_groups[id_].append(counter)
                    counter += 1
                # now do last group
                group_sel = images_in_dataset >= group_starts[-1]
                image_groups.set_selected((sel.iselection().select(group_sel)), counter)
                self.group_to_datasetid_and_range[counter] = (
                    id_,
                    (group_starts[-1], max_img),
                )
                self.datasetid_to_groups[id_].append(counter)
                counter += 1

            table["group"] = image_groups

        statistics = PerImageCChalfStatistics(
            table,
            mean_unit_cell,
            space_group,
            self.params.dmin,
            self.params.dmax,
            self.params.nbins,
        )

        statistics.run()

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

        # Remove datasets based on delta cc1/2 - if using experiments & reflections
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
                    self.group_to_datasetid_and_range,
                    self.datasetid_to_groups,
                    self.results_summary,
                )
            self.reflections = [filtered_reflections]

    def read_experiments(self, experiments, reflection_table):
        """
        Get information from experiments and reflections
        """

        # Get space group and unit cell
        space_group = None
        unit_cell_list = []
        for e in experiments:
            if space_group is None:
                space_group = e.crystal.get_space_group()
            else:
                assert (
                    space_group.type().number()
                    == e.crystal.get_space_group().type().number()
                )
            unit_cell_list.append(e.crystal.get_unit_cell())

        # Require a dials scaled experiments file.
        filtered_table = filter_reflection_table(reflection_table, ["scale"])
        filtered_table["intensity"] = filtered_table["intensity.scale.value"]
        filtered_table["variance"] = filtered_table["intensity.scale.variance"]
        filtered_table["dataset"] = filtered_table["id"]
        filtered_table["image"] = (
            flex.floor(filtered_table["xyzobs.px.value"].parts()[2]).iround() + 1
        )

        return filtered_table, unit_cell_list, space_group

    def read_mtzfile(self, filename):
        """
        Read the mtz file
        """
        miller_arrays = mtz.object(file_name=filename).as_miller_arrays(
            merge_equivalents=False
        )

        # Select the desired columns
        intensities = None
        batches = None
        for array in miller_arrays:
            if array.info().labels == ["I", "SIGI"]:
                intensities = array
            if array.info().labels == ["BATCH"]:
                batches = array
        if not intensities:
            raise KeyError("Intensities not found in mtz file, expected labels I, SIGI")
        if not batches:
            raise KeyError("Batch values not found")
        if batches.data().size() != intensities.data().size():
            raise ValueError("Batch and intensity array sizes do not match")

        # Get the unit cell and space group
        unit_cell = intensities.unit_cell()
        space_group = intensities.crystal_symmetry().space_group()

        # The reflection data
        table = flex.reflection_table()
        table["miller_index"] = intensities.indices()
        table["intensity"] = intensities.data()
        table["variance"] = intensities.sigmas() ** 2

        # Create unit cell list
        zeroed_batches = batches.data() - flex.min(batches.data())
        dataset = flex.int(table.size(), 0)
        sorted_batches = flex.sorted(zeroed_batches)
        sel_perm = flex.sort_permutation(zeroed_batches)

        batch_offset = self.params.mtz.batch_offset
        if not batch_offset:
            previous = 0
            potential_batch_offsets = flex.double()
            for i, b in enumerate(sorted_batches):
                if b - previous > 1:
                    potential_batch_offsets.append(b - previous)
                previous = b
            potential = flex.sorted(potential_batch_offsets)
            # potential is a list of low numbers (where images may not have any spots)
            # and larger numbers between batches.
            if len(potential) == 1:
                batch_offset = potential[0]
                logger.info(
                    """
Using a batch offset of %s to split datasets.
Batch offset can be specified with mtz.batch_offset=
""",
                    batch_offset,
                )
            else:
                diffs = flex.double(
                    [potential[i + 1] - p for i, p in enumerate(potential[:-1])]
                )
                i = flex.sort_permutation(diffs)[-1]
                batch_offset = int(potential[i + 1] - (0.2 * diffs[i]))
                logger.info(
                    """
Using an approximate batch offset of %s to split datasets.
Batch offset can be specified with mtz.batch_offset=
""",
                    batch_offset,
                )

        previous = 0
        dataset_no = 0
        for i, b in enumerate(sorted_batches):
            if b - previous > batch_offset - 1:
                dataset_no += 1
            dataset[i] = dataset_no
            previous = b

        table["dataset"] = flex.int(table.size(), 0)
        table["dataset"].set_selected(sel_perm, dataset)

        return table, [unit_cell], space_group

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
                if (
                    expid_to_image_groups[exp_id][-1] == id_
                    or expid_to_image_groups[exp_id][0] == id_
                ):  # is at edge of scan.
                    image_ranges_removed.append([image_range, exp_id])
                    logger.info(
                        "Removing image range %s from experiment %s",
                        image_range,
                        exp_id,
                    )
                    exclude_images.append(
                        [
                            str(exp_id)
                            + ":"
                            + str(image_range[0])
                            + ":"
                            + str(image_range[1])
                        ]
                    )
                    if expid_to_image_groups[exp_id][-1] == id_:
                        del expid_to_image_groups[exp_id][-1]
                    else:
                        del expid_to_image_groups[exp_id][0]
                    n_removed_this_cycle += 1
                else:
                    other_potential_ids_to_remove.append(id_)
            ids_to_remove = other_potential_ids_to_remove
        for id_ in other_potential_ids_to_remove:
            exp_id, image_range = image_group_to_expid_and_range[id_]
            logger.info(
                """Image range %s from experiment %s is below the cutoff, but not at the edge of a sweep.""",
                image_range,
                exp_id,
            )

        # Now remove individual batches
        if -1 in reflections["id"]:
            reflections = reflections.select(reflections["id"] != -1)
        reflection_list = reflections.split_by_experiment_id()
        reflection_list, experiments = exclude_image_ranges_for_scaling(
            reflection_list, experiments, exclude_images
        )

        # check if any image groups were all outliers and missed by the analysis
        # This catches an edge case where there is an image group full of
        # outliers, which gets filtered out before the analysis but should
        # be set as not a valid image range.
        exclude_images = []
        for i, exp in enumerate(experiments):
            if len(exp.scan.get_valid_image_ranges(exp.identifier)) > 1:
                # if any of the image ranges are not in the sets tested:
                tested = []
                for exp_id, imgrange in image_group_to_expid_and_range.values():
                    if exp_id == i:
                        tested.extend(list(range(imgrange[0], imgrange[1] + 1)))
                for imgrange in exp.scan.get_valid_image_ranges(exp.identifier):
                    if imgrange[0] not in tested:
                        exclude_images.append(
                            [str(i) + ":" + str(imgrange[0]) + ":" + str(imgrange[1])]
                        )
                        logger.info(
                            "Removing %s due to scaling outlier group.",
                            exclude_images[-1],
                        )

        if exclude_images:
            reflection_list, experiments = exclude_image_ranges_for_scaling(
                reflection_list, experiments, exclude_images
            )

        # if a whole experiment has been excluded: need to remove it here
        ids_removed = []
        for exp, refl in zip(experiments, reflection_list):
            if not exp.scan.get_valid_image_ranges(
                exp.identifier
            ):  # if all removed above
                experiments_to_delete.append(exp.identifier)
                ids_removed.append(refl.experiment_identifiers().keys()[0])
        if experiments_to_delete:
            experiments, reflection_list = select_datasets_on_identifiers(
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
                "experiment_ids_fully_removed": ids_removed,
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
        ids_removed = []
        for id_ in sorted(ids_to_remove):
            logger.info("Removing dataset %d", id_)
            datasets_to_remove.append(reflections.experiment_identifiers()[id_])
            ids_removed.append(id_)
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
                "experiment_ids_fully_removed": ids_removed,
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
            self.reflections[0].as_file(self.params.output.reflections)

            logger.info("Saving the experiments to %s", self.params.output.experiments)
            self.experiments.as_file(self.params.output.experiments)

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
    from dials.util.options import OptionParser, reflections_and_experiments_from_files

    usage = "dials.compute_delta_cchalf [options] scaled.expt scaled.refl"

    parser = OptionParser(
        usage=usage,
        phil=phil,
        epilog=help_message,
        read_experiments=True,
        read_reflections=True,
        check_format=False,
    )

    params, _ = parser.parse_args(args=args, show_diff_phil=False)

    dials.util.log.config(logfile=params.output.log)

    reflections, experiments = reflections_and_experiments_from_files(
        params.input.reflections, params.input.experiments
    )

    if len(experiments) == 0 and not params.input.mtzfile:
        parser.print_help()
        return

    script = Script(params, experiments, reflections)
    script.run()
    script.write_experiments_and_reflections()
    script.plot_data()


if __name__ == "__main__":
    with dials.util.show_mail_on_error():
        run()

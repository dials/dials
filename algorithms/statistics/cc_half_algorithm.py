"""ΔCC½ algorithm definitions"""


from __future__ import annotations

import logging
from collections import defaultdict
from math import sqrt

from jinja2 import ChoiceLoader, Environment, PackageLoader

from cctbx import uctbx
from iotbx import mtz

from dials.algorithms.scaling.scale_and_filter import (
    make_histogram_plots,
    make_per_dataset_plot,
)
from dials.algorithms.statistics.delta_cchalf import PerGroupCChalfStatistics
from dials.array_family import flex
from dials.util.exclude_images import exclude_image_ranges_for_scaling
from dials.util.filter_reflections import filter_reflection_table
from dials.util.multi_dataset_handling import select_datasets_on_identifiers

logger = logging.getLogger("dials.command_line.compute_delta_cchalf")


class CCHalfFromMTZ:

    """
    Run a cc-half algorithm using an MTZ file.
    """

    def __init__(self, params, filename):
        self.params = params

        table, unit_cell, space_group = self.read_mtzfile(
            filename, self.params.mtz.batch_offset
        )
        table["group"] = table["dataset"]

        self.algorithm = DeltaCCHalf(table, unit_cell, space_group, params)
        self.results_summary = {}

    def run(self):
        self.algorithm.run()
        # update results

        # now do the exclusion
        results = self.algorithm.results_summary
        cc_half_values = results["per_dataset_delta_cc_half_values"]
        below_cutoff = (
            flex.double(cc_half_values["delta_cc_half_values"])
            < results["dataset_removal"]["cutoff_value"]
        )
        ids_to_remove = flex.int(cc_half_values["datasets"]).select(below_cutoff)
        for id_ in sorted(ids_to_remove):
            logger.info("Dataset %d is below the cutoff", id_)
        results["dataset_removal"].update({"experiments_fully_removed": ids_to_remove})
        self.results_summary = results

    def output(self):
        self.algorithm.output_table()
        self.algorithm.output_html_report()

    @staticmethod
    def read_mtzfile(filename, batch_offset=None):
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
        table["variance"] = flex.pow2(intensities.sigmas())

        # Create unit cell list
        zeroed_batches = batches.data() - flex.min(batches.data())
        dataset = flex.int(table.size(), 0)
        sorted_batches = flex.sorted(zeroed_batches)
        sel_perm = flex.sort_permutation(zeroed_batches)

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
            elif len(potential) > 1:
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
            else:
                batch_offset = 1

        previous = 0
        dataset_no = 0
        for i, b in enumerate(sorted_batches):
            if b - previous > batch_offset - 1:
                dataset_no += 1
            dataset[i] = dataset_no
            previous = b

        table["dataset"] = flex.int(table.size(), 0)
        table["dataset"].set_selected(sel_perm, dataset)

        return table, unit_cell, space_group


class CCHalfFromDials:

    """
    Run a cc-half algorithm using dials datafiles.
    """

    def __init__(self, params, experiments, reflection_table):
        self.params = params
        self.experiments = experiments
        self.reflection_table = reflection_table
        self.filtered_reflection_table = None
        # prepare data
        self.results_summary = {}
        self.group_to_datasetid_and_range = {}
        self.datasetid_to_groups = defaultdict(list)

        table, unit_cell, space_group = self.read_experiments(
            self.experiments, self.reflection_table
        )

        if self.params.mode == "dataset":
            table["group"] = table["dataset"]
        elif self.params.mode == "image_group":
            # set up tracking of groups to expts

            image_groups = flex.int(table["dataset"].size(), 0)
            counter = 0
            for id_ in set(table["dataset"]):
                sel = table["dataset"] == id_
                expid = table.experiment_identifiers()[id_]
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
                        expid,
                        (start, group_starts[i + 1] - 1),
                    )
                    self.datasetid_to_groups[expid].append(counter)
                    counter += 1
                # now do last group
                group_sel = images_in_dataset >= group_starts[-1]
                image_groups.set_selected((sel.iselection().select(group_sel)), counter)
                self.group_to_datasetid_and_range[counter] = (
                    expid,
                    (group_starts[-1], max_img),
                )
                self.datasetid_to_groups[expid].append(counter)
                counter += 1

            table["group"] = image_groups

        self.algorithm = DeltaCCHalf(table, unit_cell, space_group, params)

    def run(self):
        """Run the ΔCC½ algorithm and then exclude data as appropriate"""
        self.algorithm.run()

        # now do the exclusion
        results = self.algorithm.results_summary
        cc_half_values = results["per_dataset_delta_cc_half_values"]
        below_cutoff = (
            flex.double(cc_half_values["delta_cc_half_values"])
            < results["dataset_removal"]["cutoff_value"]
        )
        ids_to_remove = flex.int(cc_half_values["datasets"]).select(below_cutoff)

        if self.params.mode == "dataset":
            filtered_reflections = self.remove_datasets_below_cutoff(
                self.experiments, self.reflection_table, ids_to_remove, results
            )
        elif self.params.mode == "image_group":
            filtered_reflections = self.remove_image_ranges_below_cutoff(
                self.experiments,
                self.reflection_table,
                ids_to_remove,
                self.group_to_datasetid_and_range,
                self.datasetid_to_groups,
                results,
            )
        self.filtered_reflection_table = filtered_reflections
        self.results_summary = results

    def output(self):
        """Save the output data and updated datafiles."""
        self.algorithm.output_table()
        logger.info(
            "Saving %d reflections to %s",
            self.filtered_reflection_table.size(),
            self.params.output.reflections,
        )
        self.filtered_reflection_table.as_file(self.params.output.reflections)

        logger.info("Saving the experiments to %s", self.params.output.experiments)
        self.experiments.as_file(self.params.output.experiments)
        self.algorithm.output_html_report()

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
        n_valid_reflections = reflections.get_flags(reflections.flags.scaled).count(
            True
        )

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
            output_reflections.flags.scaled
        ).count(True)
        results_summary["dataset_removal"].update(
            {
                "experiments_fully_removed": datasets_to_remove,
                "experiment_ids_fully_removed": ids_removed,
                "n_reflections_removed": n_valid_reflections
                - n_valid_filtered_reflections,
            }
        )
        return output_reflections

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
        n_valid_reflections = reflections.get_flags(reflections.flags.scaled).count(
            True
        )
        expid_to_tableid = {
            v: k
            for k, v in zip(
                reflections.experiment_identifiers().keys(),
                reflections.experiment_identifiers().values(),
            )
        }

        experiments_to_delete = []
        exclude_images = []
        image_ranges_removed = []  # track for results summary
        n_removed_this_cycle = 1
        while n_removed_this_cycle != 0:
            other_potential_ids_to_remove = []
            n_removed_this_cycle = 0
            for id_ in sorted(ids_to_remove):
                exp_id, image_range = image_group_to_expid_and_range[id_]  # identifier
                if (
                    expid_to_image_groups[exp_id][-1] == id_
                    or expid_to_image_groups[exp_id][0] == id_
                ):  # is at edge of scan.
                    # loc = list(experiments.identifiers()).index(exp_id)
                    table_id = expid_to_tableid[exp_id]
                    image_ranges_removed.append([image_range, table_id])
                    logger.info(
                        "Removing image range %s from experiment %s",
                        image_range,
                        table_id,
                    )
                    exclude_images.append(
                        [f"{table_id}:{image_range[0]}:{image_range[1]}"]
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
            table_id = expid_to_tableid[exp_id]
            logger.info(
                """Image range %s from experiment %s is below the cutoff, but not at the edge of a sweep.""",
                image_range,
                table_id,
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
        for exp in experiments:
            # if any of the image ranges are not in the sets tested, exclude them
            tested = []
            for exp_id, imgrange in image_group_to_expid_and_range.values():
                if exp_id == exp.identifier:
                    tested.extend(list(range(imgrange[0], imgrange[1] + 1)))
            for imgrange in exp.scan.get_valid_image_ranges(exp.identifier):
                if all([j not in tested for j in range(imgrange[0], imgrange[1] + 1)]):
                    table_id = expid_to_tableid[exp.identifier]
                    exclude_images.append([f"{table_id}:{imgrange[0]}:{imgrange[1]}"])
                    logger.info(
                        "Removing %s due to scaling outlier group.", exclude_images[-1]
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
            output_reflections.flags.scaled
        ).count(True)
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
    def read_experiments(experiments, reflection_table):
        """
        Get information from experiments and reflections
        """

        # Get space group and unit cell
        space_group = experiments[0].crystal.get_space_group()
        sgno = space_group.type().number()
        unit_cell_list = []
        for e in experiments:
            assert sgno == e.crystal.get_space_group().type().number()
            unit_cell_list.append(e.crystal.get_unit_cell())

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

        # Require a dials scaled experiments file.
        filtered_table = filter_reflection_table(reflection_table, ["scale"])
        filtered_table["intensity"] = filtered_table["intensity.scale.value"]
        filtered_table["variance"] = filtered_table["intensity.scale.variance"]
        filtered_table["dataset"] = filtered_table["id"]
        filtered_table["image"] = (
            flex.floor(filtered_table["xyzobs.px.value"].parts()[2]).iround() + 1
        )

        return filtered_table, mean_unit_cell, space_group


class DeltaCCHalf:

    """
    Implementation of a ΔCC½ algorithm.
    """

    def __init__(self, reflection_table, median_unit_cell, space_group, params):
        self.reflection_table = reflection_table
        self.params = params
        self.median_unit_cell = median_unit_cell
        self.space_group = space_group
        self.delta_cchalf_i = {}
        self.results_summary = {
            "dataset_removal": {
                "mode": self.params.mode,
                "stdcutoff": self.params.stdcutoff,
            }
        }
        if len(set(reflection_table["group"])) == 1:
            if self.params.mode == "dataset":
                raise ValueError(
                    """Cannot perform delta-cc-half analysis in dataset mode with only one dataset.
For image group based delta-cc-half analysis, use the option mode=image_group"""
                )
            else:
                raise ValueError(
                    """Cannot perform delta-cc-half analysis on a single image group.
Choose a suitable option for group_size to divide the dataset into multiple groups"""
                )

    def run(self):
        """Run the delta_cc_half algorithm."""

        statistics = PerGroupCChalfStatistics(
            self.reflection_table,
            self.median_unit_cell,
            self.space_group,
            self.params.dmin,
            self.params.dmax,
            self.params.nbins,
        )

        statistics.run()

        self.delta_cchalf_i = statistics.delta_cchalf_i()
        self.results_summary["mean_cc_half"] = statistics._cchalf_mean
        # Print out the datasets in order of ΔCC½
        self.sort_deltacchalf_values(self.delta_cchalf_i, self.results_summary)

        cutoff_value = self._calculate_cutoff_value(
            self.delta_cchalf_i, self.params.stdcutoff
        )
        self.results_summary["dataset_removal"].update({"cutoff_value": cutoff_value})

    def output_table(self):
        """Write a text file with delta cchalf value"""
        if self.params.output.table:
            logger.info("Writing table to %s", self.params.output.table)
            with open(self.params.output.table, "w") as outfile:
                for dataset, cchalf in zip(
                    self.results_summary["per_dataset_delta_cc_half_values"][
                        "datasets"
                    ],
                    self.results_summary["per_dataset_delta_cc_half_values"][
                        "delta_cc_half_values"
                    ],
                ):
                    outfile.write("%d %f\n" % (dataset, cchalf))

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
            logger.info("Dataset: %d, ΔCC½: %.3f", datasets[i], 100 * val)
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
        logger.info(f"\nmean delta_cc_half: {(mean * 100):.3f}")
        logger.info(f"stddev delta_cc_half: {(sdev * 100):.3f}")
        cutoff_value = mean - stdcutoff * sdev
        logger.info(f"cutoff value: {(cutoff_value * 100):.3f} \n")
        return cutoff_value

    def output_html_report(self):
        if self.params.output.html:
            data = {"cc_half_plots": {}}
            res = {
                "delta_cc_half_values": self.results_summary[
                    "per_dataset_delta_cc_half_values"
                ]["delta_cc_half_values"],
                "mean_cc_half": self.results_summary["mean_cc_half"],
            }
            if "image_ranges_removed" in self.results_summary["dataset_removal"]:
                res["image_ranges_removed"] = self.results_summary["dataset_removal"][
                    "image_ranges_removed"
                ]
            else:
                res["removed_datasets"] = self.results_summary["dataset_removal"][
                    "experiments_fully_removed"
                ]
            data["cc_half_plots"].update(make_histogram_plots([res]))
            del data["cc_half_plots"]["mean_cc_one_half_vs_cycle"]
            data["cc_half_plots"].update(make_per_dataset_plot(self.delta_cchalf_i))

            logger.info("Writing html report to: %s", self.params.output.html)
            loader = ChoiceLoader(
                [
                    PackageLoader("dials", "templates"),
                    PackageLoader("dials", "static", encoding="utf-8"),
                ]
            )
            env = Environment(loader=loader)
            template = env.get_template("simple_report.html")
            html = template.render(
                page_title="ΔCC½ report",
                panel_title="Delta CC-Half plots",
                panel_id="cc_half_plots",
                graphs=data["cc_half_plots"],
            )
            with open(self.params.output.html, "wb") as f:
                f.write(html.encode("utf-8", "xmlcharrefreplace"))

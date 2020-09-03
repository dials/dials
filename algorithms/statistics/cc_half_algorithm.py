"""ΔCC½ algorithm definitions"""

from __future__ import absolute_import, division, print_function

import logging
from collections import defaultdict
from iotbx import mtz
from cctbx import uctbx

from dials.array_family import flex
from dials.util.exclude_images import exclude_image_ranges_for_scaling
from dials.util.multi_dataset_handling import select_datasets_on_identifiers
from dials.util.filter_reflections import filter_reflection_table
from dials.algorithms.statistics.delta_cchalf import PerGroupCChalfStatistics
from dials.algorithms.scaling.scale_and_filter import (
    make_histogram_plots,
    make_per_dataset_plot,
)
from jinja2 import Environment, ChoiceLoader, PackageLoader

logger = logging.getLogger("dials.command_line.compute_delta_cchalf")


class CCHalfFromMTZ(object):

    """
    Run a cc-half algorithm using an MTZ file.
    """

    def __init__(self, params, filename):
        self.params = params

        table, unit_cell, space_group = self.read_mtzfile(
            filename, self.params.mtz.batch_offset
        )
        table["group"] = table["dataset"]

        if self.params.deltacchalf_cutoff is not None:
            cutoff = self.params.deltacchalf_cutoff
            cutoff_method = (
                "fisher" if self.params.fisher_transformation else "deltacchalf"
            )
        else:
            cutoff = self.params.stdcutoff
            cutoff_method = "normalised"
        self.statistics = PerGroupCChalfStatistics(
            table,
            unit_cell,
            space_group,
            cutoff=cutoff,
            cutoff_method=cutoff_method,
            d_min=params.d_min,
            d_max=params.d_max,
            n_bins=params.nbins,
        )

        # now do the exclusion
        self.datasets_removed = self.statistics.exclude_groups
        for id_ in sorted(self.datasets_removed):
            logger.info("Dataset %d is below the cutoff", id_)

    def output(self):
        if self.params.output.html:
            delta_cchalf_report(self, self.params.output.html)

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


class CCHalfFromDials(object):

    """
    Run a cc-half algorithm using dials datafiles.
    """

    def __init__(self, params, experiments, reflection_table):
        self.params = params
        self.experiments = experiments
        self.reflection_table = reflection_table
        self.filtered_reflection_table = None
        # prepare data
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
                min_img = flex.min(images_in_dataset)
                max_img = flex.max(images_in_dataset)
                n_groups = (max_img - min_img + 1) // self.params.group_size
                for n in range(n_groups):
                    start = min_img + n * self.params.group_size
                    if (n + 1) < n_groups:
                        end = min_img + (n + 1) * self.params.group_size - 1
                    else:
                        end = max_img
                    image_groups.set_selected(
                        (
                            sel.iselection().select(
                                (images_in_dataset >= start)
                                & (images_in_dataset <= end)
                            )
                        ),
                        counter,
                    )
                    self.group_to_datasetid_and_range[counter] = (
                        expid,
                        (start, end),
                    )
                    self.datasetid_to_groups[expid].append(counter)
                    counter += 1

            table["group"] = image_groups

        if self.params.deltacchalf_cutoff is not None:
            cutoff = self.params.deltacchalf_cutoff
            cutoff_method = (
                "fisher" if self.params.fisher_transformation else "deltacchalf"
            )
        else:
            cutoff = self.params.stdcutoff
            cutoff_method = "normalised"
        self.statistics = PerGroupCChalfStatistics(
            table,
            unit_cell,
            space_group,
            cutoff=cutoff,
            cutoff_method=cutoff_method,
            d_min=params.d_min,
            d_max=params.d_max,
            n_bins=params.nbins,
        )

        # now do the exclusion
        if self.params.mode == "dataset":
            filtered_reflections = self.remove_datasets_below_cutoff()
        elif self.params.mode == "image_group":
            filtered_reflections = self.remove_image_ranges_below_cutoff()
        self.filtered_reflection_table = filtered_reflections

    def output(self):
        """Save the output data and updated datafiles."""
        logger.info(self.statistics)
        logger.info(
            "Saving %d reflections to %s",
            self.filtered_reflection_table.size(),
            self.params.output.reflections,
        )
        self.filtered_reflection_table.as_file(self.params.output.reflections)

        logger.info("Saving the experiments to %s", self.params.output.experiments)
        self.experiments.as_file(self.params.output.experiments)

        if self.params.output.html:
            delta_cchalf_report(self, self.params.output.html)

    def remove_datasets_below_cutoff(self):
        """Remove the datasets with ids in ids_to_remove.

        Remove from the experiments and reflections and add information to the
        results summary dict.

        Returns:
          output_reflections: The reflection table with data removed.
        """
        n_valid_reflections = self.reflection_table.get_flags(
            self.reflection_table.flags.bad_for_scaling, all=False
        ).count(False)

        self.ids_removed = self.statistics.exclude_groups
        self.datasets_removed = []
        for id_ in sorted(self.ids_removed):
            logger.info("Removing dataset %d", id_)
            self.datasets_removed.append(
                self.reflection_table.experiment_identifiers()[id_]
            )
        output_reflections = self.reflection_table.remove_on_experiment_identifiers(
            self.datasets_removed
        )
        self.experiments.remove_on_experiment_identifiers(self.datasets_removed)
        output_reflections.assert_experiment_identifiers_are_consistent(
            self.experiments
        )

        n_valid_filtered_reflections = output_reflections.get_flags(
            output_reflections.flags.bad_for_scaling, all=False
        ).count(False)
        self.n_reflections_removed = n_valid_reflections - n_valid_filtered_reflections

        return output_reflections

    def remove_image_ranges_below_cutoff(self):
        """Remove image ranges from the datasets."""
        reflections = self.reflection_table.select(self.reflection_table["id"] != -1)

        n_valid_reflections = reflections.get_flags(
            reflections.flags.bad_for_scaling, all=False
        ).count(False)
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
            for id_ in sorted(self.statistics.exclude_groups):
                exp_id, image_range = self.group_to_datasetid_and_range[
                    id_
                ]  # identifier
                if (
                    self.datasetid_to_groups[exp_id][-1] == id_
                    or self.datasetid_to_groups[exp_id][0] == id_
                ):  # is at edge of scan.
                    table_id = expid_to_tableid[exp_id]
                    image_ranges_removed.append([image_range, table_id])
                    logger.info(
                        "Removing image range %s from experiment %s",
                        image_range,
                        table_id,
                    )
                    exclude_images.append(
                        ["%s:%s:%s" % (table_id, image_range[0], image_range[1])]
                    )
                    if self.datasetid_to_groups[exp_id][-1] == id_:
                        del self.datasetid_to_groups[exp_id][-1]
                    else:
                        del self.datasetid_to_groups[exp_id][0]
                    n_removed_this_cycle += 1
                else:
                    other_potential_ids_to_remove.append(id_)
            ids_removed = other_potential_ids_to_remove
        for id_ in other_potential_ids_to_remove:
            exp_id, image_range = self.group_to_datasetid_and_range[id_]
            table_id = expid_to_tableid[exp_id]
            logger.info(
                """Image range %s from experiment %s is below the cutoff, but not at the edge of a sweep.""",
                image_range,
                table_id,
            )

        # Now remove individual batches
        reflection_list = reflections.split_by_experiment_id()
        reflection_list, experiments = exclude_image_ranges_for_scaling(
            reflection_list, self.experiments, exclude_images
        )

        # check if any image groups were all outliers and missed by the analysis
        # This catches an edge case where there is an image group full of
        # outliers, which gets filtered out before the analysis but should
        # be set as not a valid image range.
        exclude_images = []
        for exp in experiments:
            # if any of the image ranges are not in the sets tested, exclude them
            tested = []
            for exp_id, imgrange in self.group_to_datasetid_and_range.values():
                if exp_id == exp.identifier:
                    tested.extend(list(range(imgrange[0], imgrange[1] + 1)))
            for imgrange in exp.scan.get_valid_image_ranges(exp.identifier):
                if all([j not in tested for j in range(imgrange[0], imgrange[1] + 1)]):
                    table_id = expid_to_tableid[exp.identifier]
                    exclude_images.append(
                        ["%s:%s:%s" % (table_id, imgrange[0], imgrange[1])]
                    )
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
            output_reflections.flags.bad_for_scaling, all=False
        ).count(False)
        self.image_ranges_removed = image_ranges_removed
        self.datasets_removed = experiments_to_delete
        self.ids_removed = ids_removed
        self.n_reflections_removed = n_valid_reflections - n_valid_filtered_reflections
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


def delta_cchalf_report(result, filename):
    data = {"cc_half_plots": {}}
    delta_cchalf = (
        result.statistics.fisher_transformed_delta_cchalf_i
        if result.statistics.cutoff_method
        else result.statistics.delta_cchalf_i
    )
    res = {
        "delta_cc_half_values": delta_cchalf,
        "mean_cc_half": result.statistics.mean_cchalf,
    }
    if result.params.mode == "image_group":
        res["image_ranges_removed"] = result.image_ranges_removed
    else:
        res["removed_datasets"] = result.datasets_removed
    res["cutoff_value"] = result.statistics.cutoff
    data = {"cc_half_plots": make_histogram_plots([res])}
    del data["cc_half_plots"]["mean_cc_one_half_vs_cycle"]
    data["cc_half_plots"].update(
        make_per_dataset_plot(result.statistics.group_ids, delta_cchalf)
    )

    logger.info("Writing html report to: %s", result.params.output.html)
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
    with open(result.params.output.html, "wb") as f:
        f.write(html.encode("utf-8", "xmlcharrefreplace"))

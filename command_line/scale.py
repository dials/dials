#!/usr/bin/env python
# coding: utf-8
from __future__ import absolute_import, division, print_function
import time
import logging
import sys
import gc
import libtbx
from libtbx import phil
from dials.util import log, show_mail_on_error, Sorry
from dials.array_family import flex
from dials.util.options import OptionParser, flatten_reflections, flatten_experiments
from dials.util.version import dials_version
from dials.algorithms.scaling.scaling_library import (
    create_scaling_model,
    create_datastructures_for_structural_model,
    create_datastructures_for_target_mtz,
    create_auto_scaling_model,
    set_image_ranges_in_scaling_models,
    scaled_data_as_miller_array,
    determine_best_unit_cell,
    merging_stats_from_scaled_array,
)
from dials.algorithms.scaling.scaler_factory import create_scaler, MultiScalerFactory
from dials.util.multi_dataset_handling import (
    select_datasets_on_ids,
    parse_multiple_datasets,
    assign_unique_identifiers,
)
from dials.util.export_mtz import match_wavelengths, make_merged_mtz_file
from dials.algorithms.scaling.scaling_utilities import (
    save_experiments,
    save_reflections,
    log_memory_usage,
    DialsMergingStatisticsError,
)
from dials.util.exclude_images import exclude_image_ranges_for_scaling
from dials.algorithms.scaling.algorithm import (
    targeted_scaling_algorithm,
    scaling_algorithm,
)
from dials.util.observer import Subject
from dials.algorithms.scaling.observers import (
    register_default_scaling_observers,
    register_merging_stats_observers,
)
from dials.report.analysis import make_merging_statistics_summary
from dials.command_line.cosym import cosym
from dials.command_line.cosym import phil_scope as cosym_phil_scope

help_message = """
This program performs scaling on integrated datasets, which attempts to improve
the internal consistency of the reflection intensities by correcting for
various experimental effects. By default, a physical scaling model is used,
with scale, decay and absorption components. If multiple input files have been
specified, the datasets will be jointly scaled against a common target of
unique reflection intensities.

The program outputs one scaled.refl and scaled.expt file, which
contains reflection data and scale models, from one or more experiments.
The output pickle file contains intensity.scale.value, the unscaled intensity
values used to determine the scaling model, and a inverse scale factor per
reflection. These values can then be used to merge the data for downstream
structural solution. Alternatively, the scaled.expt and
scaled.refl files can be passed back to dials.scale, and further scaling will
be performed, starting from where the previous job finished.

The scaling models determined by this program can be plotted with::

  dials.plot_scaling_models scaled.refl scaled.expt

Example use cases

Regular single-sweep scaling, with no absorption correction::

  dials.scale integrated.refl integrated.expt absorption_term=False

Scaling multiple datasets, specifying scale parameter interval::

  dials.scale 1_integrated.refl 1_integrated.expt 2_integrated.refl 2_integrated.expt scale_interval=10.0

Incremental scaling (with different options per dataset)::

  dials.scale integrated.refl integrated.expt scale_interval=10.0

  dials.scale integrated_2.refl integrated_2.expt scaled.refl scaled.expt scale_interval=15.0

"""


logger = logging.getLogger("dials")
phil_scope = phil.parse(
    """
  debug = False
    .type = bool
    .help = "Output additional debugging information"
  model = physical array KB
    .type = choice
    .help = "Set scaling model to be applied to input datasets without
            an existing model. "
    .expert_level = 0
  stats_only = False
    .type = bool
    .help = "Only read input files and output merging stats."
  export_mtz_only = False
    .type = bool
    .help = "Only read input scaled input files and make mtz files if"
            "user specified unmerged_mtz, merged_mtz."
  output {
    log = dials.scale.log
      .type = str
      .help = "The log filename"
    debug.log = dials.scale.debug.log
      .type = str
      .help = "The debug log filename"
    experiments = "scaled.expt"
      .type = str
      .help = "Option to set filepath for output json."
    reflections = "scaled.refl"
      .type = str
      .help = "Option to set filepath for output pickle file of scaled
               intensities."
    html = "scaling.html"
      .type = str
      .help = "Filename for html report."
    unmerged_mtz = None
      .type = strings
      .help = "Filename to export an unmerged_mtz file using dials.export."
              "Multiple names can be given to enable separated mtz export in"
              "the case of multiple wavelengths."
    merged_mtz = None
      .type = strings
      .help = "Filename to export a merged_mtz file. Multiple names can be"
              "given to enable separated mtz export in the case of multiple"
              "wavelengths."
    crystal_name = XTAL
      .type = str
      .help = "The crystal name to be exported in the mtz file metadata"
      .expert_level = 1
    use_internal_variance = False
      .type = bool
      .help = "Option to use internal spread of the intensities when merging
              reflection groups and calculating sigI, rather than using the
              sigmas of the individual reflections."
      .expert_level = 1
    merging.nbins = 20
      .type = int
      .help = "Number of bins to use for calculating and plotting merging stats."
      .expert_level = 1
    delete_integration_shoeboxes = True
      .type = bool
      .help = "Discard integration shoebox data from scaling output, to help"
              "with memory management."
      .expert_level = 2
  }
  include scope dials.algorithms.scaling.scaling_options.phil_scope
  include scope dials.algorithms.scaling.cross_validation.cross_validate.phil_scope
  include scope dials.algorithms.scaling.scaling_refiner.scaling_refinery_phil_scope
  include scope dials.util.exclude_images.phil_scope
  include scope dials.util.multi_dataset_handling.phil_scope
""",
    process_includes=True,
)


class Script(Subject):
    """Main script to run the scaling algorithm."""

    def __init__(self, params, experiments, reflections):
        super(Script, self).__init__(events=["merging_statistics", "run_script"])
        self.scaler = None
        self.scaled_miller_array = None
        self.merging_statistics_result = None
        self.anom_merging_statistics_result = None
        self.params, self.experiments, self.reflections = self.prepare_input(
            params, experiments, reflections
        )
        self._create_model_and_scaler()
        logger.debug("Initialised scaling script object")
        log_memory_usage()

    def _create_model_and_scaler(self):
        """Create the scaling models and scaler."""
        if self.params.model in (None, libtbx.Auto):
            self.experiments = create_auto_scaling_model(
                self.params, self.experiments, self.reflections
            )
        else:
            self.experiments = create_scaling_model(
                self.params, self.experiments, self.reflections
            )
        logger.info("\nScaling models have been initialised for all experiments.")
        logger.info("%s%s%s", "\n", "=" * 80, "\n")

        self.experiments = set_image_ranges_in_scaling_models(self.experiments)

        self.scaler = create_scaler(self.params, self.experiments, self.reflections)

    @Subject.notify_event(event="run_script")
    def run(self):
        """Run the scaling script."""
        start_time = time.time()
        self.scale()
        self.remove_unwanted_datasets()
        self.scaled_miller_array = scaled_data_as_miller_array(
            self.reflections, self.experiments, anomalous_flag=False
        )
        try:
            self.calculate_merging_stats()
        except DialsMergingStatisticsError as e:
            logger.info(e)

        # All done!
        logger.info("\nTotal time taken: {:.4f}s ".format(time.time() - start_time))
        logger.info("%s%s%s", "\n", "=" * 80, "\n")

    @staticmethod
    def prepare_input(params, experiments, reflections):
        """Perform checks on the data and prepare the data for scaling."""

        if params.scaling_options.space_group:
            logger.warning(
                """
Warning: the command line option space_group has been deprecated,
please use dials.reindex, dials.symmetry or dials.cosym to
prepare the data in the correct space group.\n"""
            )

        #### First exclude any datasets, before the dataset is split into
        #### individual reflection tables and expids set.
        if (
            params.dataset_selection.exclude_datasets
            or params.dataset_selection.use_datasets
        ):
            try:
                experiments, reflections = select_datasets_on_ids(
                    experiments,
                    reflections,
                    params.dataset_selection.exclude_datasets,
                    params.dataset_selection.use_datasets,
                )
                logger.info(
                    "\nDataset unique identifiers for retained datasets are %s \n",
                    list(experiments.identifiers()),
                )
            except ValueError as e:
                raise Sorry(e)

        #### Split the reflections tables into a list of reflection tables,
        #### with one table per experiment.
        logger.info(
            "Checking for the existence of a reflection table \n"
            "containing multiple datasets \n"
        )
        reflections = parse_multiple_datasets(reflections)
        logger.info(
            "Found %s reflection tables & %s experiments in total.",
            len(reflections),
            len(experiments),
        )

        if len(experiments) != len(reflections):
            raise Sorry("Mismatched number of experiments and reflection tables found.")

        #### Assign experiment identifiers.
        try:
            experiments, reflections = assign_unique_identifiers(
                experiments, reflections
            )
        except ValueError as e:
            raise Sorry(e)
        logger.info(
            "\nDataset unique identifiers are %s \n", list(experiments.identifiers())
        )
        try:
            reflections, experiments = exclude_image_ranges_for_scaling(
                reflections, experiments, params.exclude_images
            )
        except ValueError as e:
            raise Sorry(e)

        #### Allow checking of consistent indexing, useful for
        #### targeted / incremental scaling.
        if params.scaling_options.check_consistent_indexing:
            logger.info("Running dials.cosym to check consistent indexing:\n")
            cosym_params = cosym_phil_scope.extract()
            cosym_params.nproc = params.scaling_options.nproc
            cosym_instance = cosym(experiments, reflections, cosym_params)
            cosym_instance.run()
            experiments = cosym_instance.experiments
            reflections = cosym_instance.reflections
            logger.info("Finished running dials.cosym, continuing with scaling.\n")

        #### Make sure all experiments in same space group
        sgs = [e.crystal.get_space_group().type().number() for e in experiments]
        if len(set(sgs)) > 1:
            raise Sorry(
                """The experiments have different space groups:
                space group numbers found: %s
                Please reanalyse the data so that space groups are consistent,
                (consider using dials.reindex, dials.symmetry or dials.cosym) or
                remove incompatible experiments (using the option exclude_datasets=)"""
                % ", ".join(map(str, set(sgs)))
            )
        logger.info(
            "Space group being used during scaling is %s",
            experiments[0].crystal.get_space_group().info(),
        )

        #### If doing targeted scaling, extract data and append an experiment
        #### and reflection table to the lists
        if params.scaling_options.target_model:
            logger.info("Extracting data from structural model.")
            exp, reflection_table = create_datastructures_for_structural_model(
                reflections, experiments, params.scaling_options.target_model
            )
            experiments.append(exp)
            reflections.append(reflection_table)

        elif params.scaling_options.target_mtz:
            logger.info("Extracting data from merged mtz.")
            exp, reflection_table = create_datastructures_for_target_mtz(
                experiments, params.scaling_options.target_mtz
            )
            experiments.append(exp)
            reflections.append(reflection_table)

        #### Perform any non-batch cutting of the datasets, including the target dataset
        best_unit_cell = determine_best_unit_cell(experiments)
        for reflection in reflections:
            if params.cut_data.d_min or params.cut_data.d_max:
                d = best_unit_cell.d(reflection["miller_index"])
                if params.cut_data.d_min:
                    sel = d < params.cut_data.d_min
                    reflection.set_flags(sel, reflection.flags.user_excluded_in_scaling)
                if params.cut_data.d_max:
                    sel = d > params.cut_data.d_max
                    reflection.set_flags(sel, reflection.flags.user_excluded_in_scaling)
            if params.cut_data.partiality_cutoff and "partiality" in reflection:
                reflection.set_flags(
                    reflection["partiality"] < params.cut_data.partiality_cutoff,
                    reflection.flags.user_excluded_in_scaling,
                )
        return params, experiments, reflections

    def scale(self):
        """The main scaling algorithm."""

        if self.scaler.id_ == "target":
            ### FIXME add in quick prescaling round if large scale difference?
            self.scaler.perform_scaling()

            if (
                self.params.scaling_options.only_target
                or self.params.scaling_options.target_model
                or self.params.scaling_options.target_mtz
            ):

                self.scaler = targeted_scaling_algorithm(self.scaler)
                return
            # Now pass to a multiscaler ready for next round of scaling.
            self.scaler.expand_scales_to_all_reflections()
            self.scaler = MultiScalerFactory.create_from_targetscaler(self.scaler)

        # From here onwards, scaler should only be a SingleScaler
        # or MultiScaler (not TargetScaler).
        self.scaler = scaling_algorithm(self.scaler)

    def remove_unwanted_datasets(self):
        """Remove any target model/mtz data and any datasets which were removed
        from the scaler during scaling."""
        # first remove target refl/exps
        if (
            self.params.scaling_options.target_model
            or self.params.scaling_options.target_mtz
            or self.params.scaling_options.only_target
        ):
            self.experiments = self.experiments[:-1]
            self.reflections = self.reflections[:-1]

        # remove any bad datasets:
        removed_ids = self.scaler.removed_datasets
        if removed_ids:
            logger.info("deleting removed datasets from memory: %s", removed_ids)
            expids = list(self.experiments.identifiers())
            locs_in_list = []
            for id_ in removed_ids:
                locs_in_list.append(expids.index(id_))
            self.experiments, self.reflections = select_datasets_on_ids(
                self.experiments, self.reflections, exclude_datasets=removed_ids
            )

    @staticmethod
    def stats_only(reflections, experiments, params):
        """Calculate and print merging stats."""
        best_unit_cell = params.reflection_selection.best_unit_cell
        if not params.reflection_selection.best_unit_cell:
            best_unit_cell = determine_best_unit_cell(experiments)
        scaled_miller_array = scaled_data_as_miller_array(
            reflections, experiments, best_unit_cell=best_unit_cell
        )
        try:
            res, _ = merging_stats_from_scaled_array(
                scaled_miller_array,
                params.output.merging.nbins,
                params.output.use_internal_variance,
            )
            logger.info(make_merging_statistics_summary(res))
        except DialsMergingStatisticsError as e:
            logger.info(e)

    @staticmethod
    def export_mtz_only(reflections, experiments, params):
        """Export data in mtz format."""
        assert len(reflections) == 1, "Need a combined reflection table from scaling."
        if params.output.unmerged_mtz:
            _export_unmerged_mtz(params, experiments, reflections[0])

        if params.output.merged_mtz:
            if len(params.output.merged_mtz) > 1:
                _export_multi_merged_mtz(params, experiments, reflections[0])
            else:
                scaled_array = scaled_data_as_miller_array(reflections, experiments)
                merged = scaled_array.merge_equivalents(
                    use_internal_variance=params.output.use_internal_variance
                ).array()
                merged_anom = (
                    scaled_array.as_anomalous_array()
                    .merge_equivalents(
                        use_internal_variance=params.output.use_internal_variance
                    )
                    .array()
                )
                mtz_file = make_merged_mtz_file(merged, merged_anom)
                logger.info(
                    "\nSaving output to a merged mtz file to %s.\n",
                    params.output.merged_mtz[0],
                )
                mtz_file.write(params.output.merged_mtz[0])

    @Subject.notify_event(event="merging_statistics")
    def calculate_merging_stats(self):
        try:
            self.merging_statistics_result, self.anom_merging_statistics_result = merging_stats_from_scaled_array(
                self.scaled_miller_array,
                self.params.output.merging.nbins,
                self.params.output.use_internal_variance,
            )
        except DialsMergingStatisticsError as e:
            logger.info(e)

    def delete_datastructures(self):
        """Delete the data in the scaling datastructures to save RAM before
        combinining datasets for output."""
        del self.scaler
        for experiment in self.experiments:
            for component in experiment.scaling_model.components.keys():
                experiment.scaling_model.components[component] = []
        gc.collect()

    def export(self):
        """Save the experiments json and scaled pickle file."""
        logger.info("%s%s%s", "\n", "=" * 80, "\n")

        save_experiments(self.experiments, self.params.output.experiments)

        # Now create a joint reflection table. Delete all other data before
        # joining reflection tables - just need experiments for mtz export
        # and a reflection table.
        self.delete_datastructures()

        joint_table = flex.reflection_table()
        for i in range(len(self.reflections)):
            joint_table.extend(self.reflections[i])
            # del reflection_table
            self.reflections[i] = 0
            gc.collect()

        # remove reflections with neg sigma
        sel = joint_table["inverse_scale_factor"] <= 0.0
        good_sel = ~joint_table.get_flags(joint_table.flags.bad_for_scaling, all=False)
        n_neg = (good_sel & sel).count(True)
        if n_neg > 0:
            logger.warning(
                """
Warning: %s non-excluded reflections were assigned negative scale factors
during scaling. These will be set as outliers in the reflection table. It
may be best to rerun scaling from this point for an improved model.""",
                n_neg,
            )
            joint_table.set_flags(sel, joint_table.flags.outlier_in_scaling)

        save_reflections(joint_table, self.params.output.reflections)

        if self.params.output.unmerged_mtz:
            _export_unmerged_mtz(self.params, self.experiments, joint_table)

        if self.params.output.merged_mtz:
            if len(self.params.output.merged_mtz) > 1:
                _export_multi_merged_mtz(self.params, self.experiments, joint_table)
            else:
                anomalous_scaled = scaled_data_as_miller_array(
                    [joint_table], self.experiments, anomalous_flag=True
                )
                merged = self.scaled_miller_array.merge_equivalents(
                    use_internal_variance=self.params.output.use_internal_variance
                ).array()
                merged_anom = anomalous_scaled.merge_equivalents(
                    use_internal_variance=self.params.output.use_internal_variance
                ).array()
                mtz_file = make_merged_mtz_file(merged, merged_anom)
                logger.info(
                    "\nSaving output to a merged mtz file to %s.\n",
                    self.params.output.merged_mtz[0],
                )
                mtz_file.write(self.params.output.merged_mtz[0])


def _export_multi_merged_mtz(params, experiments, reflection_table):
    from dxtbx.model import ExperimentList

    wavelengths = match_wavelengths(experiments)
    assert len(params.output.merged_mtz) == len(list(wavelengths.keys()))
    for filename, wavelength in zip(params.output.merged_mtz, list(wavelengths.keys())):
        exps = ExperimentList()
        ids = []
        for i, exp in enumerate(experiments):
            if i in wavelengths[wavelength]:
                exps.append(exp)
                ids.append(exp.identifier)
        refls = reflection_table.select_on_experiment_identifiers(ids)
        scaled_array = scaled_data_as_miller_array([refls], exps)
        merged = scaled_array.merge_equivalents(
            use_internal_variance=params.output.use_internal_variance
        ).array()
        merged_anom = (
            scaled_array.as_anomalous_array()
            .merge_equivalents(
                use_internal_variance=params.output.use_internal_variance
            )
            .array()
        )
        mtz_file = make_merged_mtz_file(merged, merged_anom)
        logger.info("\nSaving output to a merged mtz file to %s.\n", filename)
        mtz_file.write(filename)


def _export_unmerged_mtz(params, experiments, reflection_table):
    """Export data to unmerged_mtz format (as single file or split by wavelength)."""
    from dials.command_line.export import MTZExporter
    from dials.command_line.export import phil_scope as export_phil_scope

    export_params = export_phil_scope.extract()

    export_params.intensity = ["scale"]
    export_params.mtz.partiality_threshold = params.cut_data.partiality_cutoff
    export_params.mtz.crystal_name = params.output.crystal_name
    if params.cut_data.d_min:
        export_params.mtz.d_min = params.cut_data.d_min
    if len(params.output.unmerged_mtz) > 1:
        from dxtbx.model import ExperimentList

        wavelengths = match_wavelengths(experiments)
        assert len(params.output.unmerged_mtz) == len(list(wavelengths.keys()))
        for filename, wavelength in zip(params.output.unmerged_mtz, list(wavelengths.keys())):
            export_params.mtz.hklout = filename
            logger.info("\nSaving output to an unmerged mtz file to %s.", filename)
            exps = ExperimentList()
            ids = []
            for i, exp in enumerate(experiments):
                if i in wavelengths[wavelength]:
                    exps.append(exp)
                    ids.append(exp.identifier)
            exporter = MTZExporter(
                export_params,
                exps,
                [reflection_table.select_on_experiment_identifiers(ids)],
            )
            exporter.export()
    else:
        logger.info(
            "\nSaving output to an unmerged mtz file to %s.",
            params.output.unmerged_mtz[0],
        )
        export_params.mtz.hklout = params.output.unmerged_mtz[0]
        exporter = MTZExporter(export_params, experiments, [reflection_table])
        exporter.export()


def run_scaling(params, experiments, reflections):
    """Run scaling algorithms; stats only, cross validation or standard."""
    if params.stats_only:
        Script.stats_only(reflections, experiments, params)
        sys.exit()

    if params.export_mtz_only:
        Script.export_mtz_only(reflections, experiments, params)
        sys.exit()

    if params.output.delete_integration_shoeboxes:
        for r in reflections:
            del r["shoebox"]

    if params.cross_validation.cross_validation_mode:
        from dials.algorithms.scaling.cross_validation.cross_validate import (
            cross_validate,
        )
        from dials.algorithms.scaling.cross_validation.crossvalidator import (
            DialsScaleCrossValidator,
        )

        cross_validator = DialsScaleCrossValidator(experiments, reflections)
        try:
            cross_validate(params, cross_validator)
        except ValueError as e:
            raise Sorry(e)

        logger.info(
            "Cross validation analysis does not produce scaling output files, rather\n"
            "it gives insight into the dataset. Choose an appropriate parameterisation\n"
            "and rerun scaling without cross_validation_mode.\n"
        )

    else:
        script = Script(params, experiments, reflections)
        # Register the observers at the highest level
        if params.output.html:
            register_default_scaling_observers(script)
        else:
            register_merging_stats_observers(script)
        script.run()
        script.export()


def run(args=None):
    """Run the scaling from the command-line."""
    usage = """Usage: dials.scale integrated.refl integrated.expt
[integrated.refl(2) integrated.expt(2) ....] [options]"""

    parser = OptionParser(
        usage=usage,
        read_experiments=True,
        read_reflections=True,
        phil=phil_scope,
        check_format=False,
        epilog=help_message,
    )
    params, _ = parser.parse_args(args=args, show_diff_phil=False)

    if not params.input.experiments or not params.input.reflections:
        parser.print_help()
        sys.exit()

    reflections = flatten_reflections(params.input.reflections)
    experiments = flatten_experiments(params.input.experiments)

    log.config(verbosity=1, info=params.output.log, debug=params.output.debug.log)
    logger.info(dials_version())

    diff_phil = parser.diff_phil.as_str()
    if diff_phil != "":
        logger.info("The following parameters have been modified:\n")
        logger.info(diff_phil)

    run_scaling(params, experiments, reflections)


if __name__ == "__main__":
    with show_mail_on_error():
        run()

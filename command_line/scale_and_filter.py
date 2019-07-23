"""
This program runs several rounds of scaling and filtering for multi-crystal
analysis, also producing summary plots.

This program runs consecutive cycles of dials.scale and dials.compute_delta_cchalf,
until one of several termination criteria are met - these are max_cycles,
max_percent_removed and min_completeness. Filtering is based on the calculation
of delta-cc-half values, that is the cc-half of the dataset when a subset of
images are removed. Any images with a highly negative cc-half are contributing
poorly to the overall merging statistics and are removed and the dataset is
then rescaled before further filtering analysis.

This program only acts to run the underlying programs and collect results in a
convenient manner - equivalent results should be obtained by running the
individual programs. All standard command-line options for these programs can
be given, alongside the termination criteria and output options for this program.
"""

from __future__ import absolute_import, division, print_function

import logging
import json
import libtbx.phil
import dials.util
from dials.command_line.scale import Script as ScalingScript
from dials.command_line.compute_delta_cchalf import Script as FilterScript
from dials.command_line.compute_delta_cchalf import phil_scope as deltacc_phil_scope
from dials.util.exclude_images import get_valid_image_ranges
from dials.util.version import dials_version
from dials.util.observer import Subject
from dials.algorithms.scaling.observers import (
    register_merging_stats_observers,
    register_scale_and_filter_observers,
)
from dials.algorithms.scaling.scale_and_filter import AnalysisResults, log_cycle_results


logger = logging.getLogger("dials")

phil_scope = libtbx.phil.parse(
    """
include scope dials.command_line.scale.phil_scope
output {
    analysis_results = "analysis_results.json"
      .type = str
      .help = "Option to set filepath for output json of analysis results."
}
""",
    process_includes=True,
)


class ScaleAndFilter(object):

    """Class to encapsulate a scaling and filtering algorithm."""

    def __init__(self, params, scaling_script, filtering_script):
        """Provide script classes and params that do scaling and filtering.

        Separate on input to (in theory) allow customisation for different
        filtering/scaling scripts."""
        super(ScaleAndFilter, self).__init__(events=["run_scale_and_filter"])
        self.scaling_script = scaling_script
        self.filtering_script = filtering_script
        self.params = params
        filter_params = deltacc_phil_scope.extract()
        filter_params.mode = params.filtering.deltacchalf.mode
        filter_params.group_size = params.filtering.deltacchalf.group_size
        filter_params.stdcutoff = params.filtering.deltacchalf.stdcutoff
        self.filtering_params = filter_params
        self.filtering_results = None

    @Subject.notify_event(event="run_scale_and_filter")
    def scale_and_filter(self, experiments, reflections):
        """Write the behaviour of the program as functions and classes outside run()"""

        results = AnalysisResults()

        for counter in range(1, self.params.filtering.deltacchalf.max_cycles + 1):
            # First run scaling, followed by filtering
            scaling_script = self.scaling_script(self.params, experiments, reflections)
            register_merging_stats_observers(scaling_script)
            scaling_script.run()
            # Log initial expids here, need to do after dataset selection in scaling
            # but before any filtering
            if counter == 1:
                results.initial_expids_and_image_ranges = [
                    (exp.identifier, exp.scan.get_image_range()) if exp.scan else None
                    for exp in experiments
                ]
            filter_script = self.filtering_script(
                self.filtering_params,
                scaling_script.experiments,
                scaling_script.reflections,
            )
            filter_script.run()

            # Reset dataset inclusion/exclusion to avoid errors for repeated scaling.
            experiments = filter_script.experiments
            reflections = filter_script.reflections
            self.params.dataset_selection.use_datasets = None
            self.params.dataset_selection.exclude_datasets = None

            valid_image_ranges = get_valid_image_ranges(experiments)
            results.expids_and_image_ranges = [
                (exp.identifier, valid_image_ranges[i]) if exp.scan else None
                for i, exp in enumerate(experiments)
            ]

            # Log results from scaling and filtering
            results = log_cycle_results(results, scaling_script, filter_script)
            logger.info(
                "Cycle %s of filtering, n_reflections removed this cycle: %s",
                counter,
                results.get_last_cycle_results()["n_removed"],
            )

            # Test termination conditions
            latest_results = results.get_last_cycle_results()
            if latest_results["n_removed"] == 0:
                logger.info(
                    "Finishing scale and filtering as no data removed in this cycle."
                )
                scaling_script.export()
                results.finish(termination_reason="no_more_removed")
                break

            if (
                latest_results["cumul_percent_removed"]
                > self.params.filtering.deltacchalf.max_percent_removed
            ):
                logger.info(
                    "Finishing scale and filtering as have now removed more than the limit."
                )
                results = self._run_final_scale(
                    self.params, experiments, reflections, results
                )
                results.finish(termination_reason="max_percent_removed")
                break

            if self.params.filtering.deltacchalf.min_completeness:
                if (
                    latest_results["merging_stats"]["completeness"]
                    < self.params.filtering.deltacchalf.min_completeness
                ):
                    logger.info(
                        "Finishing scale and filtering as completeness now below cutoff."
                    )
                    results = self._run_final_scale(
                        self.params, experiments, reflections, results
                    )
                    results.finish(termination_reason="below_completeness_limit")
                    break

            if counter == self.params.filtering.deltacchalf.max_cycles:
                logger.info("Finishing as reached max number of cycles.")
                results = self._run_final_scale(
                    self.params, experiments, reflections, results
                )
                results.finish(termination_reason="max_cycles")
                break

        # Print summary of results
        logger.info("\nSummary of data removed:")
        for i, res in enumerate(results.get_cycle_results()):
            logger.info("Cycle number: %s", i + 1)
            if "image_ranges_removed" in res:
                if res["image_ranges_removed"]:
                    logger.info(
                        "  Removed image ranges: \n    %s",
                        "\n    ".join(
                            str(t[0]) + ", dataset " + str(t[1])
                            for t in res["image_ranges_removed"]
                        ),
                    )
            else:
                if res["removed_datasets"]:
                    logger.info("  Removed datasets: %s", res["removed_datasets"])
            logger.info(
                "  cumulative %% of reflections removed: %s",
                res["cumul_percent_removed"],
            )

        self.filtering_results = results
        return results

    def _run_final_scale(self, params, experiments, reflections, results):
        scaling_script = self.scaling_script(params, experiments, reflections)
        register_merging_stats_observers(scaling_script)
        scaling_script.run()
        scaling_script.export()
        results.add_final_stats(scaling_script.merging_statistics_result)
        return results


def run(args=None, phil=phil_scope):
    """Run the scale and filter script."""
    import dials.util.log
    from dials.util.options import OptionParser
    from dials.util.options import flatten_reflections
    from dials.util.options import flatten_experiments

    usage = "dials.scale_and_filter [options] integrated.expt integrated.refl"

    parser = OptionParser(
        usage=usage,
        phil=phil,
        read_reflections=True,
        read_experiments=True,
        check_format=False,
        epilog=__doc__,
    )
    params, _ = parser.parse_args(args=args, show_diff_phil=False)

    dials.util.log.config(info=params.output.log, debug=params.output.debug.log)
    logger.info(dials_version())
    diff_phil = parser.diff_phil.as_str()
    if diff_phil:
        logger.info("The following parameters have been modified:\n%s", diff_phil)

    experiments = flatten_experiments(params.input.experiments)
    reflections = flatten_reflections(params.input.reflections)

    scale_and_filter = ScaleAndFilter(params, ScalingScript, FilterScript)
    register_scale_and_filter_observers(scale_and_filter)
    analysis_results = scale_and_filter.scale_and_filter(experiments, reflections)

    with open(params.output.analysis_results, "w") as f:
        json.dump(analysis_results.to_dict(), f, indent=2)


if __name__ == "__main__":
    with dials.util.show_mail_on_error():
        run()

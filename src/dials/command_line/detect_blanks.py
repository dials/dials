from __future__ import annotations

import json
import logging
import sys

import libtbx.phil

from dials.util import detect_blanks, show_mail_handle_errors

logger = logging.getLogger("dials.detect_blanks")

phil_scope = libtbx.phil.parse(
    """\
phi_step = 2
  .type = float(value_min=0)
  .help = "Width of bins in degrees."
counts_fractional_loss = 0.1
  .type = float(value_min=0, value_max=1)
  .help = "Fractional loss (relative to the bin with the most counts) after "
          "which a bin is flagged as potentially containing blank images."
misigma_fractional_loss = 0.1
  .type = float(value_min=0, value_max=1)
  .help = "Fractional loss (relative to the bin with the highest misigma) after "
          "which a bin is flagged as potentially containing blank images."
output {
  json = blanks.json
    .type = path
  plot = False
    .type = bool
}
""",
    process_includes=True,
)


help_message = """\
"""


@show_mail_handle_errors()
def run(args=None):
    from dials.util import log
    from dials.util.options import (
        ArgumentParser,
        reflections_and_experiments_from_files,
    )

    usage = "dials.detect_blanks [options] models.expt observations.refl"

    parser = ArgumentParser(
        usage=usage,
        phil=phil_scope,
        read_experiments=True,
        read_reflections=True,
        check_format=False,
        epilog=help_message,
    )

    params, options = parser.parse_args(args)
    reflections, experiments = reflections_and_experiments_from_files(
        params.input.reflections, params.input.experiments
    )

    if len(experiments) == 0 or len(reflections) == 0:
        parser.print_help()
        exit(0)

    # Configure the logging
    log.config(logfile="dials.detect_blanks.log")

    # Log the diff phil
    diff_phil = parser.diff_phil.as_str()
    if diff_phil != "":
        logger.info("The following parameters have been modified:\n")
        logger.info(diff_phil)

    reflections = reflections[0]

    imagesets = experiments.imagesets()

    if any(experiment.is_still() for experiment in experiments):
        sys.exit("dials.detect_blanks can only be used with rotation data")

    assert len(imagesets) == 1
    imageset = imagesets[0]
    scan = imageset.get_scan()

    integrated_sel = reflections.get_flags(reflections.flags.integrated)
    indexed_sel = reflections.get_flags(reflections.flags.indexed)
    centroid_outlier_sel = reflections.get_flags(reflections.flags.centroid_outlier)
    strong_sel = reflections.get_flags(reflections.flags.strong)
    indexed_sel &= ~centroid_outlier_sel

    logger.info(f"Analysis of {strong_sel.count(True)} strong reflections:")
    strong_results = detect_blanks.blank_counts_analysis(
        reflections.select(strong_sel),
        scan,
        phi_step=params.phi_step,
        fractional_loss=params.counts_fractional_loss,
    )
    for blank_start, blank_end in strong_results["blank_regions"]:
        logger.info(f"Potential blank images: {blank_start + 1} -> {blank_end}")

    indexed_results = None
    if indexed_sel.count(True) > 0:
        logger.info(f"Analysis of {indexed_sel.count(True)} indexed reflections:")
        indexed_results = detect_blanks.blank_counts_analysis(
            reflections.select(indexed_sel),
            scan,
            phi_step=params.phi_step,
            fractional_loss=params.counts_fractional_loss,
        )
        for blank_start, blank_end in indexed_results["blank_regions"]:
            logger.info(f"Potential blank images: {blank_start + 1} -> {blank_end}")

    integrated_results = None
    if integrated_sel.count(True) > 0:
        logger.info(f"Analysis of {integrated_sel.count(True)} integrated reflections:")
        integrated_results = detect_blanks.blank_integrated_analysis(
            reflections.select(integrated_sel),
            scan,
            phi_step=params.phi_step,
            fractional_loss=params.misigma_fractional_loss,
        )
        for blank_start, blank_end in integrated_results["blank_regions"]:
            logger.info(f"Potential blank images: {blank_start + 1} -> {blank_end}")

    d = {
        "strong": strong_results,
        "indexed": indexed_results,
        "integrated": integrated_results,
    }

    if params.output.json is not None:
        with open(params.output.json, "w") as fh:
            json.dump(d, fh)

    if params.output.plot:
        from matplotlib import pyplot

        plots = [(strong_results, "-")]
        if indexed_results:
            plots.append((indexed_results, "--"))
        if integrated_results:
            plots.append((integrated_results, ":"))

        for results, linestyle in plots:
            xs = results["data"][0]["x"]
            ys = results["data"][0]["y"]
            xmax = max(xs)
            ymax = max(ys)
            xs = [x / xmax for x in xs]
            ys = [y / ymax for y in ys]
            blanks = results["data"][0]["blank"]
            pyplot.plot(xs, ys, color="blue", linestyle=linestyle)
            pyplot.plot(
                *zip(*[(x, y) for x, y, blank in zip(xs, ys, blanks) if blank]),
                color="red",
                linestyle=linestyle,
            )
        pyplot.ylim(0)
        pyplot.show()
        pyplot.clf()


if __name__ == "__main__":
    run()

# DIALS_ENABLE_COMMAND_LINE_COMPLETION


from __future__ import annotations

import json
import logging
import pathlib
import time
from collections import Counter
from itertools import accumulate

import iotbx.phil

from dials.util import log, show_mail_handle_errors
from dials.util.options import ArgumentParser, reflections_and_experiments_from_files
from dials.util.version import dials_version

logger = logging.getLogger("dials.commandline.countrate")

help_message = """
dials.countrate: Analyse pixel intensities from reflection "shoebox" data and recommend a suitable transmission for future experiments.

dials.countrate takes the "shoeboxes" (raw pixel data) from a reflection table (e.g. strong.refl from dials.find_spots) and calculates a histogram of pixel intensities.
A maximum recommended transmission is then calculated to keep a specified percentile of pixel intensities below a target percentage of the detector trusted range. The idea
of this is to try and keep most of the pixel intensities in the linear range of the detector and limit reliance on the detector countrate correction.

Usage examples:
    dials.countrate imported.expt strong.refl
"""

phil_scope = iotbx.phil.parse(
    """
input {
    target_countrate_pct = 10.0
        .type = float
        .help = "Target percentage of the detector trusted range to scale the reference percentile reflection to"
    ref_percentile = 99.9
        .type = float
        .help = "Percentile value used in transmission recommendation. dials.countrate will calculate the transmission required to scale this percentile reflection intensity to the desired target_countrate_pct"
}

output {
    log = "dials.countrate.log"
        .type = str
        .help = "Log file for processing"
    results = "countrate.json"
        .type = str
        .help = "JSON file containing histogram of pixel intensities, detector trusted range and recommended transmission"
}
"""
)


def generate_histogram_from_shoeboxes(shoeboxes) -> dict[int, int]:
    """Process the spotfinding results and perform additional analysis."""
    logger.info("Processing spotfinding results...")

    counter: Counter[int] = Counter()
    for shoebox in shoeboxes:
        counter.update(
            int(pixel_intensity)
            for pixel_intensity in shoebox.data.as_numpy_array().ravel()
        )
    sorted_counter = sorted(counter.items())
    histogram: dict[int, int] = dict(sorted_counter)

    return histogram


def save_results_to_json(
    histogram: dict[int, int],
    max_trusted_value: int,
    results_path: pathlib.Path,
    recommended_transmission: float,
):
    logger.info(f"Saving results to {str(results_path)}\n")
    with open(results_path, "w") as f:
        json.dump(
            {
                "counts": histogram,
                "trusted_range_upper_limit": max_trusted_value,
                "recommended_transmission": recommended_transmission,
            },
            f,
            indent=2,
        )


def get_percentile_index(num_pixels: list[int], percentile: float) -> int:
    threshold = sum(num_pixels) * percentile

    for index, cum_sum in enumerate(accumulate(num_pixels)):
        if cum_sum >= threshold:
            break

    return index


@show_mail_handle_errors()
def run(args: list[str] | None = None, phil=phil_scope):
    """Main entry point for the CLI program."""
    start_time = time.time()

    # Parse command line arguments
    parser = ArgumentParser(
        usage="xia2.countrate [options] [param.phil]",
        phil=phil_scope,
        read_reflections=True,
        read_experiments=True,
        check_format=False,
        epilog=help_message,
    )

    params, options = parser.parse_args(args=args, show_diff_phil=False)

    # Setup logging
    log.config(logfile=params.output.log, verbosity=options.verbose)
    logger.info(dials_version())

    # Show parsed parameters
    diff_phil = parser.diff_phil.as_str()
    if diff_phil:
        logger.info("The following parameters have been modified:\n")
        logger.info(diff_phil)

    reflections, experiments = reflections_and_experiments_from_files(
        params.input.reflections, params.input.experiments
    )
    # Run the processing pipeline
    logger.info("=" * 60)
    logger.info("Running dials.countrate")
    logger.info("=" * 60)

    if len(experiments) != 1 or len(reflections) != 1:
        logger.error("Expected exactly one experiment and one reflection table")
        parser.print_help()
        return

    experiment = experiments[0]
    shoeboxes = reflections[0]["shoebox"]

    transmission = experiment.beam.get_transmission()

    histogram = generate_histogram_from_shoeboxes(shoeboxes)
    detector = experiment.detector.to_dict()
    max_trusted_counts = detector["panels"][0]["trusted_range"][1]

    max_pixel_count = max(histogram.keys())
    max_pixel_percent_of_trusted_range = max_pixel_count * 100 / max_trusted_counts

    num_pixels = list(histogram.values())
    pixel_intensities = list(histogram.keys())
    total_pixels = sum(num_pixels)

    percentiles = [99.999, 99.99, 99.9, 99.0, 90.0]
    percentile_trusted_range_pct: list[float] = []
    for percentile in percentiles:
        percentile_idx = get_percentile_index(num_pixels, percentile / 100)
        percentile_counts = pixel_intensities[percentile_idx]
        percentile_trusted_range_pct.append(
            percentile_counts * 100 / max_trusted_counts
        )

    # Calculate transmission needed to scale reference percentile reflection to target countrate
    target_countrate_pct = params.input.target_countrate_pct
    target_counts = max_trusted_counts * (target_countrate_pct / 100)
    ref_percentile = params.input.ref_percentile
    ref_percentile_idx = get_percentile_index(num_pixels, ref_percentile / 100)
    ref_percentile_counts = pixel_intensities[ref_percentile_idx]
    scale_factor = target_counts / ref_percentile_counts
    recommended_transmission = min(transmission * scale_factor, 1.0)

    # Log summary
    logger.info("=" * 60)
    logger.info("Processing Summary:")
    logger.info(f"Found {total_pixels} pixels from {len(shoeboxes)} reflections\n")
    logger.info(f"Experiment transmission = {transmission * 100:.2f}%")
    logger.info(
        f"Max pixel recorded at {max_pixel_percent_of_trusted_range:.2f}% of detector trusted range\n"
    )
    for percentile_index, percentile in enumerate(percentiles):
        logger.info(
            f"{percentile}% of pixels <= {percentile_trusted_range_pct[percentile_index]:.2f}% of detector trusted range\n"
        )
    logger.info(
        f"Recommended max transmission of {recommended_transmission * 100:.2f}% to keep {ref_percentile}% of pixel intensities below {target_countrate_pct}% of detector trusted range\n"
    )
    logger.info("=" * 60)

    save_results_to_json(
        histogram, max_trusted_counts, params.output.results, recommended_transmission
    )

    duration = time.time() - start_time

    logger.info(
        f"Processing took {time.strftime('%Hh %Mm %Ss', time.gmtime(duration))}"
    )


if __name__ == "__main__":
    run()

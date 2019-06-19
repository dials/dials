from __future__ import absolute_import, division, print_function

import logging
import math

import libtbx.phil
from libtbx.math_utils import iceil
from dials.util import Sorry
from scitbx.array_family import flex

logger = logging.getLogger("dials.command_line.detect_blanks")

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


def blank_counts_analysis(reflections, scan, phi_step, fractional_loss):
    if not len(reflections):
        raise Sorry("Input contains no reflections")

    xyz_px = reflections["xyzobs.px.value"]
    x_px, y_px, z_px = xyz_px.parts()
    phi = scan.get_angle_from_array_index(z_px)

    osc = scan.get_oscillation()[1]
    n_images_per_step = iceil(phi_step / osc)
    phi_step = n_images_per_step * osc

    array_range = scan.get_array_range()
    phi_min = scan.get_angle_from_array_index(array_range[0])
    phi_max = scan.get_angle_from_array_index(array_range[1])
    assert phi_min <= flex.min(phi)
    assert phi_max >= flex.max(phi)
    n_steps = int(round((phi_max - phi_min) / phi_step))
    hist = flex.histogram(
        z_px, data_min=array_range[0], data_max=array_range[1], n_slots=n_steps
    )
    logger.debug("Histogram:")
    logger.debug(hist.as_str())

    counts = hist.slots()
    fractional_counts = counts.as_double() / flex.max(counts)

    potential_blank_sel = fractional_counts <= fractional_loss

    xmin, xmax = zip(
        *[
            (slot_info.low_cutoff, slot_info.high_cutoff)
            for slot_info in hist.slot_infos()
        ]
    )

    d = {
        "data": [
            {
                "x": list(hist.slot_centers()),
                "y": list(hist.slots()),
                "xlow": xmin,
                "xhigh": xmax,
                "blank": list(potential_blank_sel),
                "type": "bar",
                "name": "blank_counts_analysis",
            }
        ],
        "layout": {
            "xaxis": {"title": "z observed (images)"},
            "yaxis": {"title": "Number of reflections"},
            "bargap": 0,
        },
    }

    blank_regions = blank_regions_from_sel(d["data"][0])
    d["blank_regions"] = blank_regions

    return d


def blank_integrated_analysis(reflections, scan, phi_step, fractional_loss):
    prf_sel = reflections.get_flags(reflections.flags.integrated_prf)
    if prf_sel.count(True) > 0:
        reflections = reflections.select(prf_sel)
        intensities = reflections["intensity.prf.value"]
        variances = reflections["intensity.prf.variance"]
    else:
        sum_sel = reflections.get_flags(reflections.flags.integrated_sum)
        reflections = reflections.select(sum_sel)
        intensities = reflections["intensity.sum.value"]
        variances = reflections["intensity.sum.variance"]

    i_sigi = intensities / flex.sqrt(variances)

    xyz_px = reflections["xyzobs.px.value"]
    x_px, y_px, z_px = xyz_px.parts()
    phi = scan.get_angle_from_array_index(z_px)

    osc = scan.get_oscillation()[1]
    n_images_per_step = iceil(phi_step / osc)
    phi_step = n_images_per_step * osc

    array_range = scan.get_array_range()
    phi_min = flex.min(phi)
    phi_max = flex.max(phi)
    n_steps = int(round((phi_max - phi_min) / phi_step))
    hist = flex.histogram(
        z_px, data_min=array_range[0], data_max=array_range[1], n_slots=n_steps
    )
    logger.debug("Histogram:")
    logger.debug(hist.as_str())

    mean_i_sigi = flex.double()
    for i, slot_info in enumerate(hist.slot_infos()):
        sel = (z_px >= slot_info.low_cutoff) & (z_px < slot_info.high_cutoff)
        if sel.count(True) == 0:
            mean_i_sigi.append(0)
        else:
            mean_i_sigi.append(flex.mean(i_sigi.select(sel)))

    potential_blank_sel = mean_i_sigi <= (fractional_loss * flex.max(mean_i_sigi))

    xmin, xmax = zip(
        *[
            (slot_info.low_cutoff, slot_info.high_cutoff)
            for slot_info in hist.slot_infos()
        ]
    )

    d = {
        "data": [
            {
                "x": list(hist.slot_centers()),
                "y": list(mean_i_sigi),
                "xlow": xmin,
                "xhigh": xmax,
                "blank": list(potential_blank_sel),
                "type": "bar",
                "name": "blank_counts_analysis",
            }
        ],
        "layout": {
            "xaxis": {"title": "z observed (images)"},
            "yaxis": {"title": "Number of reflections"},
            "bargap": 0,
        },
    }

    blank_regions = blank_regions_from_sel(d["data"][0])
    d["blank_regions"] = blank_regions

    return d


def blank_regions_from_sel(d):
    blank_sel = d["blank"]
    xlow = d["xlow"]
    xhigh = d["xhigh"]

    blank_regions = []
    n = len(blank_sel)

    for i in range(len(blank_sel)):
        if blank_sel[i]:
            if i == 0 or not blank_sel[i - 1]:
                blank_start = math.floor(xlow[i])
            blank_end = math.ceil(xhigh[i])
        if (not blank_sel[i] and i > 0 and blank_sel[i - 1]) or (
            i == (n - 1) and blank_sel[i - 1]
        ):
            blank_regions.append((blank_start, blank_end))

    return blank_regions


def run(args):

    from dials.util.options import OptionParser
    from dials.util.options import flatten_experiments
    from dials.util.options import flatten_reflections
    from dials.util import log
    import libtbx.load_env

    usage = "dials.detect_blanks [options] models.expt observations.refl"

    parser = OptionParser(
        usage=usage,
        phil=phil_scope,
        read_experiments=True,
        read_reflections=True,
        check_format=False,
        epilog=help_message,
    )

    params, options = parser.parse_args()
    experiments = flatten_experiments(params.input.experiments)
    reflections = flatten_reflections(params.input.reflections)

    if len(experiments) == 0 or len(reflections) == 0:
        parser.print_help()
        exit(0)

    # Configure the logging
    log.config(info="dials.detect_blanks.log")

    # Log the diff phil
    diff_phil = parser.diff_phil.as_str()
    if diff_phil != "":
        logger.info("The following parameters have been modified:\n")
        logger.info(diff_phil)

    reflections = reflections[0]

    imagesets = experiments.imagesets()

    assert len(imagesets) == 1
    imageset = imagesets[0]
    scan = imageset.get_scan()

    integrated_sel = reflections.get_flags(reflections.flags.integrated)
    indexed_sel = reflections.get_flags(reflections.flags.indexed)
    centroid_outlier_sel = reflections.get_flags(reflections.flags.centroid_outlier)
    strong_sel = reflections.get_flags(reflections.flags.strong)
    indexed_sel &= ~centroid_outlier_sel

    logger.info("Analysis of %i strong reflections:" % strong_sel.count(True))
    strong_results = blank_counts_analysis(
        reflections.select(strong_sel),
        scan,
        phi_step=params.phi_step,
        fractional_loss=params.counts_fractional_loss,
    )
    for blank_start, blank_end in strong_results["blank_regions"]:
        logger.info("Potential blank images: %i -> %i" % (blank_start + 1, blank_end))

    indexed_results = None
    if indexed_sel.count(True) > 0:
        logger.info("Analysis of %i indexed reflections:" % indexed_sel.count(True))
        indexed_results = blank_counts_analysis(
            reflections.select(indexed_sel),
            scan,
            phi_step=params.phi_step,
            fractional_loss=params.counts_fractional_loss,
        )
        for blank_start, blank_end in indexed_results["blank_regions"]:
            logger.info(
                "Potential blank images: %i -> %i" % (blank_start + 1, blank_end)
            )

    integrated_results = None
    if integrated_sel.count(True) > 0:
        logger.info(
            "Analysis of %i integrated reflections:" % integrated_sel.count(True)
        )
        integrated_results = blank_integrated_analysis(
            reflections.select(integrated_sel),
            scan,
            phi_step=params.phi_step,
            fractional_loss=params.misigma_fractional_loss,
        )
        for blank_start, blank_end in integrated_results["blank_regions"]:
            logger.info(
                "Potential blank images: %i -> %i" % (blank_start + 1, blank_end)
            )

    d = {
        "strong": strong_results,
        "indexed": indexed_results,
        "integrated": integrated_results,
    }

    if params.output.json is not None:
        import json

        with open(params.output.json, "wb") as f:
            json.dump(d, f)

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
                linestyle=linestyle
            )
        pyplot.ylim(0)
        pyplot.show()
        pyplot.clf()


if __name__ == "__main__":
    import sys

    run(sys.argv[1:])

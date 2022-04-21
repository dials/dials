from __future__ import annotations

import logging
import math

from libtbx.math_utils import iceil

from dials.array_family import flex

logger = logging.getLogger(__name__)


def blank_counts_analysis(reflections, scan, phi_step, fractional_loss):
    if not len(reflections):
        raise ValueError("Input contains no reflections")

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
    n_steps = max(int(round((phi_max - phi_min) / phi_step)), 1)
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

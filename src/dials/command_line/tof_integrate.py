# LIBTBX_SET_DISPATCHER_NAME dials.tof_integrate
from __future__ import annotations

import logging
import multiprocessing

import numpy as np

import cctbx.array_family.flex
import libtbx
from dxtbx import flumpy

import dials.util.log
from dials.algorithms.integration.integrator import (
    generate_phil_scope as integrator_phil_scope,
)
from dials.algorithms.integration.report import IntegrationReport
from dials.algorithms.shoebox import MaskCode
from dials.algorithms.spot_prediction import TOFReflectionPredictor
from dials.array_family import flex
from dials.command_line.integrate import process_reference
from dials.extensions.simple_background_ext import SimpleBackgroundExt
from dials.util.options import ArgumentParser, reflections_and_experiments_from_files
from dials.util.phil import parse
from dials.util.version import dials_version
from dials_algorithms_tof_integration_ext import (
    TOFProfile1DParams,
    TOFProfile3DParams,
    integrate_reflection_table,
    tof_calculate_ellipse_shoebox_mask,
    tof_calculate_seed_skewness_shoebox_mask,
)
from dials_tof_scaling_ext import (
    TOFAbsorptionParams,
    TOFIncidentSpectrumParams,
    tof_extract_shoeboxes_to_reflection_table,
)

logger = logging.getLogger("dials.command_line.simple_integrate")

phil_scope = parse(
    """
output {
experiments = 'integrated.expt'
    .type = str
    .help = "The experiments output filename"
reflections = 'integrated.refl'
    .type = str
    .help = "The integrated output filename"
phil = 'dials.simple_integrate.phil'
    .type = str
    .help = "The output phil file"
log = 'tof_integrate.log'
    .type = str
    .help = "The log filename"
}
integration_type = *observed calculated
    .type = choice
    .help = "observed integrates only reflections observed during spot finding"
            "calculated integrates all reflections out to calculated.dmin"

background_model = constant2d constant3d linear2d *linear3d
    .type = choice
    .help = Model used to calculate the background of each reflection

calculated{
    dmin = none
    .type = float (value_min=0.5)
    .help = "The resolution spots are integrated to when using integration_type.calculated"

}
method = *summation profile1d profile3d
    .type = choice
    .help = "Integration method: "
            "summation: shoebox summation"
            "profile1d: https://doi.org/10.1038/srep36628 "
            "profile3d: https://doi.org/10.1016/j.nima.2016.12.026"

mask = *ellipse seed_skewness
    .type = choice
    .help = "Foreground/background mask method: "
            "seed_skewness: https://doi.org/10.1107/S0021889803021939"


corrections{
    incident_run = None
        .type = str
        .help = "Path to incident run to normalize intensities (e.g. Vanadium)."
    empty_run = None
        .type = str
        .help = "Path to empty run to correct empty counts."
    lorentz = False
        .type = bool
        .help = "Apply the Lorentz correction to target spectrum."
    absorption{
        incident_spectrum{
            sample_number_density = 0.0722
                .type = float
                .help = "Sample number density for incident run."
                        "Default is Vanadium used at SXD"
            sample_radius = 0.03
                .type = float
                .help = "Sample radius incident run."
                        "Default is Vanadium used at SXD"
            scattering_x_section = 5.158
                .type = float
                .help = "Sample scattering cross section used for incident run."
                        "Default is Vanadium used at SXD"
            absorption_x_section = 4.4883
                .type = float
                .help = "Sample absorption cross section for incident run."
                        "Default is Vanadium used at SXD"
        }
        target_spectrum{
            sample_number_density = None
                .type = float
                .help = "Sample number density for target run."
            sample_radius = None
                .type = float
                .help = "Sample radius target run."
            scattering_x_section = None
                .type = float
                .help = "Sample scattering cross section used for target run."
            absorption_x_section = None
                .type = float
                .help = "Sample absorption cross section for target run."
        }
    }
}
profile1d{
    init_alpha = 0.03
        .type = float
        .help = "Initial alpha value before optimization"
    init_beta = 0.03
        .type = float
        .help = "Initial beta value before optimization"
    min_alpha = 0.02
        .type = float
        .help = "Min alpha value for optimization"
    max_alpha = 1.0
        .type = float
        .help = "Max alpha value for optimization"
    min_beta = 0.0
        .type = float
        .help = "Min beta value for optimization"
    max_beta = 1.0
        .type = float
        .help = "Max beta value for optimization"
    n_restarts = 8
        .type = int(value_min=0)
        .help = "If fit fails, number of additional attempts with perturbed params"

}
profile3d{
    init_alpha = 1.0
        .type = float
        .help = "Initial alpha value before optimization"
    init_beta = 0.1
        .type = float
        .help = "Initial beta value before optimization"
    min_alpha = 0.5
        .type = float
        .help = "Min alpha value for optimization"
    max_alpha = 20.
        .type = float
        .help = "Max alpha value for optimization"
    min_beta = 1e-7
        .type = float
        .help = "Min beta value for optimization"
    max_beta = 5.0
        .type = float
        .help = "Max beta value for optimization"
    n_restarts = 8
        .type = int(value_min=0)
        .help = "If fit fails, number of additional attempts with perturbed params"
    gradient_method = *forward_difference central_difference
        .type = choice
        .help = "Method used to calculate gradients"
}

mp{
    nproc = auto
        .type = int(value_min=1)
        .help = "Number of processors to use during parallelized steps."
        "If set to auto DIALS will choose automatically."
}
bbox_tof_padding = 2
    .type = int
    .help = "Additional ToF frames added to calculated bounding boxes"
bbox_xy_padding = 1
    .type = int
    .help = "Additional pixels added to calculated bounding boxes"
keep_shoeboxes = False
    .type = bool
    .help = "Retain shoeboxes in output reflection table"
"""
)

"""
Kabsch 2010 refers to
Kabsch W., Integration, scaling, space-group assignment and
post-refinment, Acta Crystallographica Section D, 2010, D66, 133-144
Usage:
$ dials.tof_integrate.py refined.expt refined.refl
"""

phil_scope.adopt_scope(integrator_phil_scope())


def get_corrections_data(experiments, params):
    corrections = {}

    experiment_cls = experiments[0].imageset.get_format_class()
    incident_proton_charge = None
    empty_proton_charge = None

    if applying_incident_and_empty_runs(params):
        incident_fmt_class = experiment_cls.get_instance(
            params.corrections.incident_run
        )
        incident_data = experiment_cls(params.corrections.incident_run).get_imageset(
            params.corrections.incident_run
        )
        corrections["incident_data"] = incident_data

        empty_data = experiment_cls(params.corrections.empty_run).get_imageset(
            params.corrections.empty_run
        )
        corrections["empty_data"] = empty_data
        empty_fmt_class = experiment_cls.get_instance(params.corrections.empty_run)
        incident_proton_charge = incident_fmt_class.get_proton_charge()
        empty_proton_charge = empty_fmt_class.get_proton_charge()

        corrections["incident_proton_charge"] = incident_proton_charge
        corrections["empty_proton_charge"] = empty_proton_charge

        if applying_spherical_absorption_correction(params):
            corrections["absorption_params"] = TOFAbsorptionParams(
                params.corrections.absorption.target_spectrum.sample_radius,
                params.corrections.absorption.target_spectrum.scattering_x_section,
                params.corrections.absorption.target_spectrum.absorption_x_section,
                params.corrections.absorption.target_spectrum.sample_number_density,
                params.corrections.absorption.incident_spectrum.sample_radius,
                params.corrections.absorption.incident_spectrum.scattering_x_section,
                params.corrections.absorption.incident_spectrum.absorption_x_section,
                params.corrections.absorption.incident_spectrum.sample_number_density,
            )

    return corrections


def calculate_shoebox_masks(experiment, reflections, method):
    if method == "seed_skewness":
        logger.info("    Calculating seed skewness foreground/background mask")
        tof_calculate_seed_skewness_shoebox_mask(reflections, experiment, 1e-7, 10)
    else:
        logger.info("    Calculating ellipse foreground/background mask")
        tof_calculate_ellipse_shoebox_mask(reflections, experiment)

    return reflections


def integrate_reflection_table_for_experiment(
    expt, expt_reflections, expt_data, params, **kwargs
):
    apply_lorentz = params.corrections.lorentz
    profile1d_params = None
    profile3d_params = None
    incident_params = None
    absorption_params = None

    logger.info(f"    Integrating using {params.method}")

    if params.method == "profile1d":
        alpha = params.profile1d.init_alpha
        beta = params.profile1d.init_beta
        A = 1.0
        min_alpha = params.profile1d.min_alpha
        max_alpha = params.profile1d.max_alpha
        min_beta = params.profile1d.min_beta
        max_beta = params.profile1d.max_beta
        n_restarts = params.profile1d.n_restarts
        profile1d_params = TOFProfile1DParams(
            A, alpha, min_alpha, max_alpha, beta, min_beta, max_beta, n_restarts, True
        )
    elif params.method == "profile3d":
        alpha = params.profile3d.init_alpha
        beta = params.profile3d.init_beta
        min_alpha = params.profile3d.min_alpha
        max_alpha = params.profile3d.max_alpha
        min_beta = params.profile3d.min_beta
        max_beta = params.profile3d.max_beta
        n_restarts = params.profile3d.n_restarts
        use_central_diff = params.profile3d.gradient_method == "central_difference"
        profile3d_params = TOFProfile3DParams(
            alpha,
            min_alpha,
            max_alpha,
            beta,
            min_beta,
            max_beta,
            n_restarts,
            True,
            use_central_diff,
        )

    if apply_lorentz:
        logger.info("    Adding Lorentz correction")

    if "incident_data" in kwargs:
        incident_params = TOFIncidentSpectrumParams(
            kwargs["incident_data"],
            kwargs["empty_data"],
            kwargs["expt_proton_charge"],
            kwargs["incident_proton_charge"],
            kwargs["empty_proton_charge"],
        )
        logger.info("    Adding incident spectrum correction")

        if "absorption_params" in kwargs:
            logger.info("    Adding absorption correction")
            absorption_params = kwargs["absorption_params"]

    integrate_reflection_table(
        expt_reflections,
        expt,
        expt_data,
        incident_params,
        absorption_params,
        apply_lorentz,
        params.mp.nproc,
        profile1d_params,
        profile3d_params,
    )

    return expt_reflections


def remove_overlapping_reflections(reflections):
    overlaps = reflections.find_overlaps()
    overlap_sel = flex.bool(len(reflections), False)
    for item in overlaps.edges():
        i0 = overlaps.source(item)
        i1 = overlaps.target(item)
        overlap_sel[i0] = True
        overlap_sel[i1] = True
    logger.info("Rejecting %i overlapping bounding boxes", overlap_sel.count(True))
    reflections = reflections.select(~overlap_sel)
    return reflections


def compute_partiality(bbox, image_size):
    intersect_x0 = max(image_size[0], bbox[0])
    intersect_y0 = max(image_size[2], bbox[2])
    intersect_x1 = min(image_size[1], bbox[1])
    intersect_y1 = min(image_size[3], bbox[3])

    if intersect_x0 >= intersect_x1 or intersect_y0 >= intersect_y1:
        return 0

    intersection_area = (intersect_x1 - intersect_x0) * (intersect_y1 - intersect_y0)
    square_area = (bbox[1] - bbox[0]) * (bbox[3] - bbox[2])

    return intersection_area / square_area


def update_bounding_box(bbox, centroid, new_centroid, padding, image_size):
    from copy import deepcopy
    from math import ceil, floor

    diff_centroid = (
        new_centroid[0] - centroid[0],
        new_centroid[1] - centroid[1],
        new_centroid[2] - centroid[2],
    )

    updated_bbox = list(deepcopy(bbox))

    updated_bbox[0] += diff_centroid[0] - padding[0]
    updated_bbox[1] += diff_centroid[0] + padding[0]
    updated_bbox[2] += diff_centroid[1] - padding[1]
    updated_bbox[3] += diff_centroid[1] + padding[1]
    updated_bbox[4] += diff_centroid[2] - padding[2]
    updated_bbox[5] += diff_centroid[2] + padding[2]

    partiality = compute_partiality(updated_bbox, image_size)

    updated_bbox[0] = max(floor(updated_bbox[0]), image_size[0])
    updated_bbox[1] = min(ceil(updated_bbox[1]), image_size[1])
    updated_bbox[2] = max(floor(updated_bbox[2]), image_size[2])
    updated_bbox[3] = min(ceil(updated_bbox[3]), image_size[3])
    updated_bbox[4] = max(floor(updated_bbox[4]), image_size[4])
    updated_bbox[5] = min(ceil(updated_bbox[5]), image_size[5])

    return tuple(updated_bbox), partiality


def split_reflections(reflections, n, by_panel=False):
    if by_panel:
        for i in range(max(reflections["panel"]) + 1):
            sel = reflections["panel"] == i
            yield reflections.select(sel)
    else:
        d, r = divmod(len(reflections), n)
        for i in range(n):
            si = (d + 1) * (i if i < r else r) + d * (0 if i < r else i - r)
            yield reflections[si : si + (d + 1 if i < r else d)]


def join_reflections(list_of_reflections):
    reflections = list_of_reflections[0]
    for i in range(1, len(list_of_reflections)):
        reflections.extend(list_of_reflections[i])
    return reflections


def get_predicted_observed_reflections(params, experiments, reflections):
    min_s0_idx = min(
        range(len(reflections["wavelength"])), key=reflections["wavelength"].__getitem__
    )
    min_s0 = reflections["s0"][min_s0_idx]
    dmin = None
    for experiment in experiments:
        expt_dmin = experiment.detector.get_max_resolution(min_s0)
        if dmin is None or expt_dmin < dmin:
            dmin = expt_dmin

    logger.info(f"Calculated dmin from observed reflections: {round(dmin, 3)}")
    predicted_reflections = None
    miller_indices = reflections["miller_index"]
    for idx, experiment in enumerate(experiments):
        predictor = TOFReflectionPredictor(experiment, dmin)
        if predicted_reflections is None:
            predicted_reflections = predictor.indices_for_ub(miller_indices)
            predicted_reflections["id"] = cctbx.array_family.flex.int(
                len(predicted_reflections), idx
            )
            predicted_reflections["imageset_id"] = cctbx.array_family.flex.int(
                len(predicted_reflections), idx
            )
        else:
            r = predictor.indices_for_ub(miller_indices)
            r["id"] = cctbx.array_family.flex.int(len(r), idx)
            r["imageset_id"] = cctbx.array_family.flex.int(len(r), idx)
            predicted_reflections.extend(r)
    predicted_reflections["s0"] = predicted_reflections["s0_cal"]
    predicted_reflections.calculate_entering_flags(experiments)

    for i in range(len(experiments)):
        predicted_reflections.experiment_identifiers()[i] = experiments[i].identifier

    # Updates flags to set which reflections to use in generating reference profiles
    matched, reflections, unmatched = predicted_reflections.match_with_reference(
        reflections
    )
    # sel = predicted_reflections.get_flags(predicted_reflections.flags.reference_spot)
    predicted_reflections = predicted_reflections.select(matched)
    if "idx" in reflections:
        predicted_reflections["idx"] = reflections["idx"]

    predicted_reflections["xyzobs.px.value"] = reflections["xyzobs.px.value"]

    tof_padding = params.bbox_tof_padding
    xy_padding = params.bbox_xy_padding
    image_size = experiments[0].detector[0].get_image_size()
    tof_size = len(experiments[0].scan.get_property("time_of_flight"))
    bboxes = flex.int6(len(predicted_reflections))
    partiality = flex.double(len(predicted_reflections))
    for i in range(len(predicted_reflections)):
        bboxes[i], partiality[i] = update_bounding_box(
            reflections["bbox"][i],
            reflections["xyzobs.px.value"][i],
            predicted_reflections["xyzcal.px"][i],
            (int(xy_padding), int(xy_padding), int(tof_padding)),
            (0, image_size[0], 0, image_size[1], 0, tof_size),
        )
    predicted_reflections["bbox"] = bboxes
    predicted_reflections["partiality"] = partiality

    predicted_reflections.compute_d(experiments)
    predicted_reflections["shoebox"] = flex.shoebox(
        predicted_reflections["panel"],
        predicted_reflections["bbox"],
        allocate=False,
        flatten=False,
    )

    predicted_reflections.map_centroids_to_reciprocal_space(
        experiments, calculated=True
    )
    return predicted_reflections


def get_predicted_calculated_reflections(params, experiments, reflections):
    dmin = params.calculated.dmin
    assert dmin is not None, (
        "Integrating calculated reflections but calculated.dmin has not been set"
    )

    predicted_reflections = None
    for idx, experiment in enumerate(experiments):
        predictor = TOFReflectionPredictor(experiment, dmin)
        if predicted_reflections is None:
            predicted_reflections = predictor.for_ub(experiment.crystal.get_A())
            predicted_reflections["id"] = cctbx.array_family.flex.int(
                len(predicted_reflections), idx
            )
            predicted_reflections["imageset_id"] = cctbx.array_family.flex.int(
                len(predicted_reflections), idx
            )
        else:
            r = predictor.for_ub(experiment.crystal.get_A())
            r["id"] = cctbx.array_family.flex.int(len(r), idx)
            r["imageset_id"] = cctbx.array_family.flex.int(len(r), idx)
            predicted_reflections.extend(r)
    predicted_reflections["s0"] = predicted_reflections["s0_cal"]
    predicted_reflections.calculate_entering_flags(experiments)

    for i in range(len(experiments)):
        predicted_reflections.experiment_identifiers()[i] = experiments[i].identifier

    logger.info(f"Predicted {len(predicted_reflections)} using dmin {dmin}")

    used_in_ref = reflections.get_flags(reflections.flags.used_in_refinement)
    model_reflections = reflections.select(used_in_ref)
    logger.info(
        f"Using {len(model_reflections)} observed reflections to calculate bounding boxes"
    )

    # Get the closest observed reflection for each calculated reflection
    # Use the observed bounding box as a basis for the calculated reflection

    tof_padding = params.bbox_tof_padding
    xy_padding = params.bbox_xy_padding
    image_size = experiments[0].detector[0].get_image_size()
    tof_size = len(experiments[0].scan.get_property("time_of_flight"))
    predicted_reflections["bbox"] = flex.int6(len(predicted_reflections))
    predicted_reflections["partiality"] = flex.double(len(predicted_reflections))

    # If no observed reflections on a given panel, use the average bbox size
    x0, x1, y0, y1, z0, z1 = reflections["bbox"].parts()
    num_reflections = len(reflections)
    avg_bbox_size = (
        (sum(x1 - x0) / num_reflections) * 0.5,
        (sum(y1 - y0) / num_reflections) * 0.5,
        (sum(z1 - z0) / num_reflections) * 0.5,
    )

    for panel_idx in range(len(experiments[0].detector)):
        for expt_idx in range(len(experiments)):
            p_sel = (predicted_reflections["id"] == expt_idx) & (
                predicted_reflections["panel"] == panel_idx
            )
            expt_p_reflections = predicted_reflections.select(p_sel)
            bboxes = flex.int6(len(expt_p_reflections))
            partiality = flex.double(len(expt_p_reflections))

            m_sel = (model_reflections["id"] == expt_idx) & (
                model_reflections["panel"] == panel_idx
            )
            expt_m_reflections = model_reflections.select(m_sel)

            for i in range(len(expt_p_reflections)):
                if len(expt_m_reflections) == 0:
                    c = expt_p_reflections["xyzcal.px"][i]
                    bbox = (
                        c[0] - avg_bbox_size[0],
                        c[0] + avg_bbox_size[0],
                        c[1] - avg_bbox_size[1],
                        c[1] + avg_bbox_size[1],
                        c[2] - avg_bbox_size[2],
                        c[2] + avg_bbox_size[2],
                    )

                    bboxes[i], partiality[i] = update_bounding_box(
                        bbox,
                        expt_p_reflections["xyzcal.px"][i],
                        expt_p_reflections["xyzcal.px"][i],
                        (int(xy_padding), int(xy_padding), int(tof_padding)),
                        (0, image_size[0], 0, image_size[1], 0, tof_size),
                    )

                else:
                    distances = (
                        expt_m_reflections["xyzobs.px.value"]
                        - expt_p_reflections["xyzcal.px"][i]
                    ).norms()
                    min_distance_idx = np.argmin(flumpy.to_numpy(distances))

                    bboxes[i], partiality[i] = update_bounding_box(
                        expt_m_reflections["bbox"][min_distance_idx],
                        expt_m_reflections["xyzobs.px.value"][min_distance_idx],
                        expt_p_reflections["xyzcal.px"][i],
                        (int(xy_padding), int(xy_padding), int(tof_padding)),
                        (0, image_size[0], 0, image_size[1], 0, tof_size),
                    )
            predicted_reflections["bbox"].set_selected(p_sel, bboxes)
            predicted_reflections["partiality"].set_selected(p_sel, partiality)

    _, _, _, _, z1, z2 = predicted_reflections["bbox"].parts()
    sel = z2 > z1
    predicted_reflections = predicted_reflections.select(sel)
    _, _, y1, y2, _, _ = predicted_reflections["bbox"].parts()
    sel = y2 > y1
    predicted_reflections = predicted_reflections.select(sel)
    x1, x2, _, _, _, _ = predicted_reflections["bbox"].parts()
    sel = x2 > x1
    predicted_reflections = predicted_reflections.select(sel)

    predicted_reflections.compute_d(experiments)
    predicted_reflections["shoebox"] = flex.shoebox(
        predicted_reflections["panel"],
        predicted_reflections["bbox"],
        allocate=False,
        flatten=False,
    )

    predicted_reflections.map_centroids_to_reciprocal_space(
        experiments, calculated=True
    )
    return predicted_reflections


def run():
    """
    Input setup
    """

    phil = phil_scope.fetch()

    usage = "usage: dials.tof_integrate.py refined.expt refined.refl"
    parser = ArgumentParser(
        usage=usage,
        phil=phil,
        epilog=__doc__,
        read_experiments=True,
        read_reflections=True,
    )

    params, options = parser.parse_args(args=None, show_diff_phil=False)

    dials.util.log.config(verbosity=options.verbose, logfile=params.output.log)
    logger.info(dials_version())

    """
    Load experiment and reflections
    """

    reflections, experiments = reflections_and_experiments_from_files(
        params.input.reflections, params.input.experiments
    )
    reflections = reflections[0]

    integrated_reflections = run_integrate(params, experiments, reflections)
    integrated_reflections.as_msgpack_file(params.output.reflections)
    experiments.as_file(params.output.experiments)


def applying_spherical_absorption_correction(params):
    all_params_present = True
    some_params_present = False
    for i in dir(params.corrections.absorption.target_spectrum):
        if i.startswith("__"):
            continue
        if getattr(params.corrections.absorption.target_spectrum, i) is not None:
            some_params_present = True
        else:
            all_params_present = False
    if some_params_present and not all_params_present:
        raise ValueError(
            "Trying to apply spherical absorption correction but some corrections are None."
        )
    return all_params_present


def applying_incident_and_empty_runs(params):
    if params.corrections.incident_run is not None:
        assert params.corrections.empty_run is not None, (
            "Incident run given without empty run."
        )
        return True
    elif params.corrections.empty_run is not None:
        assert params.corrections.incident_run is not None, (
            "Empty run given without incident run."
        )
        return True
    return False


def run_integrate(params, experiments, reflections):
    if params.mp.nproc is libtbx.Auto:
        params.mp.nproc = multiprocessing.cpu_count()

    reflections, _ = process_reference(reflections)

    """
    Predict reflections using experiment crystal
    """

    if params.integration_type == "observed":
        predicted_reflections = get_predicted_observed_reflections(
            params=params, experiments=experiments, reflections=reflections
        )
    elif params.integration_type == "calculated":
        predicted_reflections = get_predicted_calculated_reflections(
            params=params, experiments=experiments, reflections=reflections
        )

    predicted_reflections = remove_overlapping_reflections(predicted_reflections)

    corrections_data = get_corrections_data(experiments=experiments, params=params)

    experiment_cls = experiments[0].imageset.get_format_class()
    for idx, expt in enumerate(experiments):
        logger.info(f"Integrating experiment {idx}")

        expt_data = expt.imageset

        if "incident_proton_charge" in corrections_data:
            expt_proton_charge = experiment_cls.get_instance(
                expt.imageset.paths()[0], **expt.imageset.data().get_params()
            ).get_proton_charge()
            corrections_data["expt_proton_charge"] = expt_proton_charge

        sel_expt = predicted_reflections["id"] == idx
        expt_reflections = predicted_reflections.select(sel_expt)

        logger.info("    Extracting shoeboxes")
        tof_extract_shoeboxes_to_reflection_table(
            expt_reflections, expt, expt_data, False
        )

        expt_reflections = calculate_shoebox_masks(expt, expt_reflections, params.mask)
        expt_reflections.is_overloaded(experiments)
        expt_reflections.contains_invalid_pixels()

        # Filter reflections with a high fraction of masked foreground
        valid_foreground_threshold = 1.0  # DIALS default
        sboxs = expt_reflections["shoebox"]
        nvalfg = sboxs.count_mask_values(MaskCode.Valid | MaskCode.Foreground)
        nforeg = sboxs.count_mask_values(MaskCode.Foreground)
        fraction_valid = nvalfg.as_double() / nforeg.as_double()
        selection = fraction_valid < valid_foreground_threshold
        expt_reflections.set_flags(selection, expt_reflections.flags.dont_integrate)

        expt_reflections["num_pixels.valid"] = sboxs.count_mask_values(MaskCode.Valid)
        expt_reflections["num_pixels.background"] = sboxs.count_mask_values(
            MaskCode.Valid | MaskCode.Background
        )
        expt_reflections["num_pixels.background_used"] = sboxs.count_mask_values(
            MaskCode.Valid | MaskCode.Background | MaskCode.BackgroundUsed
        )
        expt_reflections["num_pixels.foreground"] = nvalfg

        logger.info(f"    Calculating background using {params.background_model}")
        params.integration.background.simple.model.algorithm = params.background_model
        background_algorithm = SimpleBackgroundExt(
            params=params, experiments=experiments
        )
        success = background_algorithm.compute_background(expt_reflections)
        expt_reflections.set_flags(
            ~success, expt_reflections.flags.failed_during_background_modelling
        )

        expt_reflections = integrate_reflection_table_for_experiment(
            expt,
            expt_reflections,
            expt_data,
            params,
            **corrections_data,
        )

        predicted_reflections.set_selected(sel_expt, expt_reflections)

    integration_report = IntegrationReport(experiments, predicted_reflections)
    logger.info("")
    logger.info(integration_report.as_str(prefix=" "))

    if not params.keep_shoeboxes:
        del predicted_reflections["shoebox"]
    return predicted_reflections


if __name__ == "__main__":
    run()

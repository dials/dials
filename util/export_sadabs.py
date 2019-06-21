from __future__ import absolute_import, division, print_function

import os
import re

import logging
from dials.util.filter_reflections import filter_reflection_table

logger = logging.getLogger(__name__)


def export_sadabs(integrated_data, experiment_list, params):
    """Export data from integrated_data corresponding to experiment_list to a
    file for input to SADABS. FIXME probably need to make a .p4p file as
    well..."""

    from dials.array_family import flex
    from scitbx import matrix
    import math

    # for the moment assume (and assert) that we will convert data from exactly
    # one lattice...

    assert len(experiment_list) == 1
    # select reflections that are assigned to an experiment (i.e. non-negative id)

    integrated_data = integrated_data.select(integrated_data["id"] >= 0)
    assert max(integrated_data["id"]) == 0

    # export for sadabs should only be for non-scaled reflections
    assert any(
        [i in integrated_data for i in ["intensity.sum.value", "intensity.prf.value"]]
    )

    integrated_data = filter_reflection_table(
        integrated_data,
        intensity_choice=params.intensity,
        partiality_threshold=params.mtz.partiality_threshold,
        combine_partials=params.mtz.combine_partials,
        min_isigi=params.mtz.min_isigi,
        filter_ice_rings=params.mtz.filter_ice_rings,
        d_min=params.mtz.d_min,
    )

    experiment = experiment_list[0]
    assert not experiment.scan is None

    # sort data before output
    nref = len(integrated_data["miller_index"])
    indices = flex.size_t_range(nref)
    perm = sorted(indices, key=lambda k: integrated_data["miller_index"][k])
    integrated_data = integrated_data.select(flex.size_t(perm))

    assert not experiment.goniometer is None

    # Warn of unhelpful SADABS behaviour for certain multi-sweep data sets
    hkl_file_root, _ = os.path.splitext(params.sadabs.hklout)
    if not params.sadabs.run or re.search("_0+$", hkl_file_root):
        logger.warning(
            "\nWarning:\n"
            "It seems SADABS rejects multi-sweep data when the first "
            "filename ends "
            "'_0', '_00', etc., with a cryptic error message:\n"
            "\t'Inconsistent 2theta values in same scan'.\n"
            "You may need to begin the numbering of your SADABS HKL files from 1, "
            "rather than 0, and ensure the SADABS run/batch number is greater than 0.\n"
        )

    axis = matrix.col(experiment.goniometer.get_rotation_axis_datum())

    beam = matrix.col(experiment.beam.get_direction())
    s0 = matrix.col(experiment.beam.get_s0())

    F = matrix.sqr(experiment.goniometer.get_fixed_rotation())
    S = matrix.sqr(experiment.goniometer.get_setting_rotation())
    unit_cell = experiment.crystal.get_unit_cell()

    if params.debug:
        m_format = "%6.3f%6.3f%6.3f\n%6.3f%6.3f%6.3f\n%6.3f%6.3f%6.3f"
        c_format = "%.2f %.2f %.2f %.2f %.2f %.2f"

        logger.info(
            "Unit cell parameters from experiment: %s"
            % (c_format % unit_cell.parameters())
        )
        logger.info(
            "Symmetry: %s" % experiment.crystal.get_space_group().type().lookup_symbol()
        )

        logger.info("Goniometer fixed matrix:\n%s" % (m_format % F.elems))
        logger.info("Goniometer setting matrix:\n%s" % (m_format % S.elems))
        logger.info("Goniometer scan axis:\n%6.3f%6.3f%6.3f" % (axis.elems))

    # detector scaling info
    assert len(experiment.detector) == 1
    panel = experiment.detector[0]
    dims = panel.get_image_size()
    pixel = panel.get_pixel_size()
    fast_axis = matrix.col(panel.get_fast_axis())
    slow_axis = matrix.col(panel.get_slow_axis())
    normal = fast_axis.cross(slow_axis)
    detector2t = s0.angle(normal, deg=True)
    origin = matrix.col(panel.get_origin())

    if params.debug:
        logger.info("Detector fast, slow axes:")
        logger.info("%6.3f%6.3f%6.3f" % (fast_axis.elems))
        logger.info("%6.3f%6.3f%6.3f" % (slow_axis.elems))
        logger.info("Detector two theta (degrees): %.2f" % detector2t)

    scl_x = 512.0 / (dims[0] * pixel[0])
    scl_y = 512.0 / (dims[1] * pixel[1])

    image_range = experiment.scan.get_image_range()

    from cctbx.array_family import flex as cflex  # implicit import # noqa: F401
    from cctbx.miller import map_to_asu_isym  # implicit import # noqa: F401

    # gather the required information for the reflection file

    nref = len(integrated_data["miller_index"])
    zdet = flex.double(integrated_data["xyzcal.px"].parts()[2])

    miller_index = integrated_data["miller_index"]

    if "intensity.sum.value" in integrated_data:
        I = integrated_data["intensity.sum.value"]
        V = integrated_data["intensity.sum.variance"]
        assert V.all_gt(0)
        sigI = flex.sqrt(V)
    else:
        I = integrated_data["intensity.prf.value"]
        V = integrated_data["intensity.prf.variance"]
        assert V.all_gt(0)
        sigI = flex.sqrt(V)

    # figure out scaling to make sure data fit into format 2F8.2 i.e. Imax < 1e5

    Imax = flex.max(I)

    if params.debug:
        logger.info("Maximum intensity in file: %8.2f" % Imax)

    if Imax > 99999.0:
        scale = 99999.0 / Imax
        I = I * scale
        sigI = sigI * scale

    phi_start, phi_range = experiment.scan.get_image_oscillation(image_range[0])

    if params.sadabs.predict:
        logger.info("Using scan static predicted spot locations")
        from dials.algorithms.spot_prediction import ScanStaticReflectionPredictor

        predictor = ScanStaticReflectionPredictor(experiment)
        UB = experiment.crystal.get_A()
        predictor.for_reflection_table(integrated_data, UB)

    if not experiment.crystal.num_scan_points:
        logger.info("No scan varying model: use static")
        static = True
    else:
        static = False

    with open(params.sadabs.hklout, "w") as fout:

        for j in range(nref):

            h, k, l = miller_index[j]

            if params.sadabs.predict:
                x_mm, y_mm, z_rad = integrated_data["xyzcal.mm"][j]
            else:
                x_mm, y_mm, z_rad = integrated_data["xyzobs.mm.value"][j]

            z0 = integrated_data["xyzcal.px"][j][2]
            istol = int(round(10000 * unit_cell.stol((h, k, l))))

            if params.sadabs.predict or static:
                # work from a scan static model & assume perfect goniometer
                # FIXME maybe should work back in the option to predict spot positions
                UB = matrix.sqr(experiment.crystal.get_A())
                phi = phi_start + z0 * phi_range
                R = axis.axis_and_angle_as_r3_rotation_matrix(phi, deg=True)
                RUB = S * R * F * UB
            else:
                # properly compute RUB for every reflection
                UB = matrix.sqr(experiment.crystal.get_A_at_scan_point(int(round(z0))))
                phi = phi_start + z0 * phi_range
                R = axis.axis_and_angle_as_r3_rotation_matrix(phi, deg=True)
                RUB = S * R * F * UB

            x = RUB * (h, k, l)
            s = (s0 + x).normalize()

            # can also compute s based on centre of mass of spot
            # s = (origin + x_mm * fast_axis + y_mm * slow_axis).normalize()

            astar = (RUB * (1, 0, 0)).normalize()
            bstar = (RUB * (0, 1, 0)).normalize()
            cstar = (RUB * (0, 0, 1)).normalize()

            ix = beam.dot(astar)
            iy = beam.dot(bstar)
            iz = beam.dot(cstar)

            dx = s.dot(astar)
            dy = s.dot(bstar)
            dz = s.dot(cstar)

            x = x_mm * scl_x
            y = y_mm * scl_y
            z = (z_rad * 180 / math.pi - phi_start) / phi_range

            fout.write(
                "%4d%4d%4d%8.2f%8.2f%4d%8.5f%8.5f%8.5f%8.5f%8.5f%8.5f"
                % (h, k, l, I[j], sigI[j], params.sadabs.run, ix, dx, iy, dy, iz, dz)
            )
            fout.write("%7.2f%7.2f%8.2f%7.2f%5d\n" % (x, y, z, detector2t, istol))

    fout.close()
    logger.info("Output %d reflections to %s" % (nref, params.sadabs.hklout))

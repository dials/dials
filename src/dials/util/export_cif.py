from __future__ import annotations

import datetime
import logging
import math

import gemmi

from cctbx import crystal, miller, uctbx
from cctbx.eltbx import tiny_pse
from cctbx.eltbx.formula import formula
from iotbx.merging_statistics import dataset_statistics
from libtbx.utils import format_float_with_standard_uncertainty

from dials.algorithms.scaling.scaling_library import determine_best_unit_cell
from dials.array_family import flex
from dials.util import Sorry
from dials.util.export_shelx import (
    parse_compound,
    single_wavelength,
    unit_cell_and_esds,
    unit_cell_volume_esd,
)
from dials.util.filter_reflections import filter_reflection_table
from dials.util.version import dials_version

logger = logging.getLogger(__name__)

DIALS_CITATION = "DIALS (Winter, G. W. et al., 2018)"


def export_cif(scaled_data, experiment_list, params):
    """Export scaled data corresponding to experiment_list to a small molecule
    (core) CIF file, suitable for structure solution in e.g. Olex2."""

    if experiment_list.all_laue() or experiment_list.all_tof():
        raise Sorry(
            "Polychromatic data cannot yet be exported as CIF. Use format=shelx instead."
        )

    # Handle requesting profile intensities (default via auto) but no column
    if "profile" in params.intensity and "intensity.prf.value" not in scaled_data:
        raise Sorry(
            "Requested profile intensity data but only summed present. Use intensity=sum."
        )

    # Check that all experiments have the same space group
    if len({x.crystal.get_space_group().make_tidy() for x in experiment_list}) != 1:
        raise ValueError("Experiments do not have a unique space group")

    # use supplied best unit cell or that determined from experiment list to
    # define d in the reflection table
    best_unit_cell = params.mtz.best_unit_cell
    if best_unit_cell is None:
        best_unit_cell = determine_best_unit_cell(experiment_list)
    else:
        logger.info("Using supplied unit cell across experiments : %s", best_unit_cell)
    scaled_data["d"] = best_unit_cell.d(scaled_data["miller_index"])

    # Clean up the reflection table with mtz defaults (as in export_shelx)
    scaled_data = filter_reflection_table(
        scaled_data,
        intensity_choice=params.intensity,
        partiality_threshold=params.mtz.partiality_threshold,
        combine_partials=params.mtz.combine_partials,
        min_isigi=params.mtz.min_isigi,
        filter_ice_rings=params.mtz.filter_ice_rings,
        d_min=params.mtz.d_min,
    )

    space_group = experiment_list[0].crystal.get_space_group()
    intensity_choice = (
        params.intensity[0] if params.intensity[0] != "profile" else "prf"
    )

    # Olex2 treats a 0 0 0 index as "ignore everything from here on", so make
    # sure that no such reflection is written
    sel = scaled_data["miller_index"] != (0, 0, 0)
    if sel.count(False):
        logger.info("Removing %d reflections with index 0 0 0", sel.count(False))
        scaled_data = scaled_data.select(sel)
    if not scaled_data.size():
        raise Sorry("No reflections left to export")

    intensities = scaled_data[f"intensity.{intensity_choice}.value"]
    variances = scaled_data[f"intensity.{intensity_choice}.variance"]
    assert variances.all_gt(0)
    sigmas = flex.sqrt(variances)

    miller_set = miller.set(
        crystal_symmetry=crystal.symmetry(
            unit_cell=best_unit_cell, space_group=space_group
        ),
        indices=scaled_data["miller_index"],
        anomalous_flag=False,
    )
    i_obs = miller_set.array(data=intensities)
    i_obs.set_observation_type_xray_intensity()
    i_obs.set_sigmas(sigmas)

    wavelength = single_wavelength(experiment_list)

    doc = gemmi.cif.Document()
    block = doc.add_new_block(params.cif.datablock)

    _write_audit(block)
    _write_chemical(block, params.cif.chemical_formula, space_group)
    _write_cell(block, experiment_list, params, scaled_data, wavelength)
    _write_symmetry(block, space_group)
    _write_experiment(block, experiment_list, wavelength)
    _write_reflection_statistics(block, i_obs, scaled_data, wavelength)
    _write_absorption(block, experiment_list, params)
    # written last, so that a user-supplied value wins over a derived one
    _write_extra(block, params.cif.extra)
    _write_reflections(block, scaled_data, intensities, sigmas, params)

    options = gemmi.cif.WriteOptions()
    options.align_pairs = 34
    options.align_loops = 16
    doc.write_file(params.cif.hklout, options)
    logger.info("Wrote %d reflections to %s", scaled_data.size(), params.cif.hklout)


def _write_audit(block):
    block.set_pair("_audit_creation_method", gemmi.cif.quote(dials_version()))
    block.set_pair("_audit_creation_date", datetime.date.today().isoformat())
    block.set_pair("_computing_data_reduction", gemmi.cif.quote(DIALS_CITATION))
    block.set_pair("_computing_cell_refinement", gemmi.cif.quote(DIALS_CITATION))


def _write_chemical(block, chemical_formula, space_group):
    """Write the chemical formula items, if a composition was supplied. A
    composition given as a bare list of elements carries no counts, so only the
    formula string can be written in that case."""

    if chemical_formula is None:
        return
    counts = parse_compound(chemical_formula)
    if not any(counts.values()):
        logger.info(
            "Composition %s has no element counts, so only _chemical_formula_sum "
            "will be written",
            chemical_formula,
        )
        block.set_pair(
            "_chemical_formula_sum", gemmi.cif.quote(" ".join(counts.keys()))
        )
        return

    block.set_pair(
        "_chemical_formula_sum",
        gemmi.cif.quote(
            str(formula(counts).sorted_as_c_h_then_by_increasing_atomic_number())
        ),
    )
    weight = sum(n * tiny_pse.table(e).weight() for e, n in counts.items())
    block.set_pair("_chemical_formula_weight", f"{weight:.2f}")
    block.set_pair("_cell_formula_units_Z", str(space_group.order_z()))


def _constrained_angles(space_group, unit_cell):
    """Return a 3-tuple of flags marking which of alpha, beta and gamma are
    fixed by the symmetry of the space group, and so should not be given a
    standard uncertainty."""

    flags = []
    for i in range(3, 6):
        params = list(unit_cell.parameters())
        params[i] += 2.0
        flags.append(not space_group.is_compatible_unit_cell(uctbx.unit_cell(params)))
    return flags


def _write_cell(block, experiment_list, params, reflections, wavelength):
    uc, uc_sd = unit_cell_and_esds(experiment_list, params.mtz.best_unit_cell)
    space_group = experiment_list[0].crystal.get_space_group()
    fixed = (False, False, False) + tuple(_constrained_angles(space_group, uc))

    names = (
        "length_a",
        "length_b",
        "length_c",
        "angle_alpha",
        "angle_beta",
        "angle_gamma",
    )
    for name, value, sd, is_fixed in zip(
        names, uc.parameters(), uc_sd or (None,) * 6, fixed
    ):
        if is_fixed or sd is None or sd <= 0:
            block.set_pair(f"_cell_{name}", f"{value:.4f}")
        else:
            block.set_pair(
                f"_cell_{name}", format_float_with_standard_uncertainty(value, sd)
            )

    volume_sd = None if uc_sd is None else unit_cell_volume_esd(experiment_list)
    if volume_sd is None:
        block.set_pair("_cell_volume", f"{uc.volume():.2f}")
    else:
        block.set_pair(
            "_cell_volume",
            format_float_with_standard_uncertainty(uc.volume(), volume_sd),
        )

    # The cell was refined against all the reflections being exported
    d = reflections["d"]
    block.set_pair("_cell_measurement_reflns_used", str(reflections.size()))
    block.set_pair(
        "_cell_measurement_theta_min", f"{_theta_from_d(max(d), wavelength):.4f}"
    )
    block.set_pair(
        "_cell_measurement_theta_max", f"{_theta_from_d(min(d), wavelength):.4f}"
    )
    block.set_pair("_cell_measurement_wavelength", f"{wavelength:.5f}")


def _theta_from_d(d, wavelength):
    """Bragg angle theta in degrees for a resolution d and wavelength."""

    return math.degrees(math.asin(min(wavelength / (2.0 * d), 1.0)))


def _write_symmetry(block, space_group):
    sg_type = space_group.type()
    loop = block.init_loop("_space_group_symop_", ["id", "operation_xyz"])
    for i, op in enumerate(space_group):
        loop.add_row([str(i + 1), gemmi.cif.quote(op.as_xyz())])
    block.set_pair(
        "_space_group_crystal_system",
        gemmi.cif.quote(space_group.crystal_system().lower()),
    )
    block.set_pair("_space_group_IT_number", str(sg_type.number()))
    block.set_pair(
        "_space_group_name_H-M_alt", gemmi.cif.quote(sg_type.lookup_symbol())
    )
    block.set_pair(
        "_space_group_name_Hall", gemmi.cif.quote(sg_type.hall_symbol().strip())
    )
    # Aliases from the superseded symmetry category, still expected by some programs
    block.set_pair(
        "_symmetry_space_group_name_H-M", gemmi.cif.quote(sg_type.lookup_symbol())
    )
    block.set_pair("_symmetry_Int_Tables_number", str(sg_type.number()))


def electron_voltage(wavelength):
    """The accelerating voltage in volts that gives electrons of this de Broglie
    wavelength in Angstroms, by relativistic inversion of

        lambda = h / sqrt(2 m0 e V (1 + e V / (2 m0 c^2)))
    """

    # CODATA 2018 values, SI units
    h = 6.62607015e-34
    m0 = 9.1093837015e-31
    e = 1.602176634e-19
    c = 299792458.0

    # a quadratic in V, of which only the positive root is physical
    a = 2 * m0 * e * e / (2 * m0 * c * c)
    b = 2 * m0 * e
    d = -((h / (wavelength * 1e-10)) ** 2)
    return (-b + math.sqrt(b * b - 4 * a * d)) / (2 * a)


def _rotation_mode(experiment_list):
    """Classify the data collection as continuous rotation or stepwise, in the
    sense of _diffrn_measurement_rotation_mode."""

    for exp in experiment_list:
        if exp.scan is None or exp.scan.get_oscillation()[1] == 0.0:
            # stills, or a scan of zero width per image
            return "stepwise"
    return "rotation"


def _write_experiment(block, experiment_list, wavelength):
    probe = {
        "x-ray": "x-ray",
        "xray": "x-ray",
        "electron": "electron",
        "neutron": "neutron",
    }.get(experiment_list[0].beam.get_probe_name(), "x-ray")
    block.set_pair("_diffrn_radiation_probe", probe)
    block.set_pair("_diffrn_radiation_wavelength", f"{wavelength:.5f}")

    if probe == "electron":
        _write_electron_diffraction(block, experiment_list, wavelength)
    elif all(exp.goniometer is not None for exp in experiment_list):
        block.set_pair("_diffrn_measurement_method", gemmi.cif.quote("\\w scans"))

    block.set_pair("_diffrn_detector", gemmi.cif.quote("area detector"))
    detector_type = experiment_list[0].detector[0].get_type()
    if detector_type:
        block.set_pair("_diffrn_detector_type", gemmi.cif.quote(detector_type))

    # The orientation matrix is only meaningful for a single experiment
    if len(experiment_list) == 1:
        block.set_pair(
            "_diffrn_orient_matrix_type", gemmi.cif.quote("DIALS convention")
        )
        A = experiment_list[0].crystal.get_A()
        for i in range(3):
            for j in range(3):
                block.set_pair(
                    f"_diffrn_orient_matrix_UB_{i + 1}{j + 1}", f"{A[3 * i + j]:.7f}"
                )


def _write_electron_diffraction(block, experiment_list, wavelength):
    """Write the electron diffraction specific items of the core dictionary.

    Only those that follow from the experimental models are written. The rest
    of the items described by the electron diffraction extension, such as the
    precession semi-angle, the illumination mode or the make of the gun, are
    not knowable from a DIALS data reduction and must be supplied by the user
    with cif.extra (or added downstream by e.g. Olex2)."""

    rotation_mode = _rotation_mode(experiment_list)

    block.set_pair("_diffrn_source", gemmi.cif.quote("electron gun"))
    block.set_pair(
        "_diffrn_source_voltage", f"{electron_voltage(wavelength) / 1e3:.0f}"
    )
    block.set_pair(
        "_diffrn_measurement_device",
        gemmi.cif.quote("transmission electron microscope"),
    )
    block.set_pair("_diffrn_measurement_rotation_mode", rotation_mode)
    method = (
        "continuous rotation 3D electron diffraction"
        if rotation_mode == "rotation"
        else "stepwise 3D electron diffraction"
    )
    block.set_pair("_diffrn_measurement_method", gemmi.cif.quote(method))


def _write_extra(block, extra):
    """Write user-supplied '_data_name=value' pairs, which override anything
    already written for the same data name."""

    for item in extra:
        name, sep, value = item.partition("=")
        name = name.strip()
        value = value.strip()
        if not sep or not name.startswith("_"):
            raise Sorry(
                f"Cannot interpret cif.extra={item}. Expected a CIF data name "
                "and value, such as cif.extra=_diffrn_source_type='LaB6 gun'"
            )
        # the value may already be quoted, in which case do not quote it again
        if len(value) > 1 and value[0] == value[-1] and value[0] in "'\"":
            value = value[1:-1]
        if value in ("?", "."):
            # the CIF special values for unknown and not applicable, which are
            # only meaningful unquoted
            block.set_pair(name, value)
        else:
            block.set_pair(name, gemmi.cif.quote(value))


def _write_reflection_statistics(block, i_obs, reflections, wavelength):
    stats = dataset_statistics(
        i_obs=i_obs,
        crystal_symmetry=i_obs.crystal_symmetry(),
        use_internal_variance=False,
        eliminate_sys_absent=False,
        assert_is_not_unique_set_under_symmetry=False,
    )
    overall = stats.overall

    block.set_pair("_diffrn_reflns_number", str(reflections.size()))
    if overall.r_merge is not None:
        block.set_pair("_diffrn_reflns_av_R_equivalents", f"{overall.r_merge:.4f}")
    sum_i = sum(abs(v) for v in i_obs.data())
    if sum_i > 0:
        block.set_pair(
            "_diffrn_reflns_av_unetI/netI", f"{sum(i_obs.sigmas()) / sum_i:.4f}"
        )

    span = miller.index_span(reflections["miller_index"])
    for name, lo, hi in zip("hkl", span.min(), span.max()):
        block.set_pair(f"_diffrn_reflns_limit_{name}_min", str(lo))
        block.set_pair(f"_diffrn_reflns_limit_{name}_max", str(hi))

    d = reflections["d"]
    theta_min = _theta_from_d(max(d), wavelength)
    theta_max = _theta_from_d(min(d), wavelength)
    block.set_pair("_diffrn_reflns_theta_min", f"{theta_min:.4f}")
    block.set_pair("_diffrn_reflns_theta_max", f"{theta_max:.4f}")
    block.set_pair("_diffrn_reflns_theta_full", f"{theta_max:.4f}")
    block.set_pair("_diffrn_measured_fraction_theta_max", f"{overall.completeness:.4f}")
    block.set_pair(
        "_diffrn_measured_fraction_theta_full", f"{overall.completeness:.4f}"
    )
    block.set_pair(
        "_diffrn_reflns_reduction_process",
        gemmi.cif.quote(_merging_summary(overall)),
    )

    merged = i_obs.merge_equivalents(use_internal_variance=False).array()
    block.set_pair("_reflns_number_total", str(merged.size()))
    block.set_pair(
        "_reflns_number_gt", str((merged.data() > 2 * merged.sigmas()).count(True))
    )
    block.set_pair("_reflns_threshold_expression", gemmi.cif.quote("I > 2\\s(I)"))
    block.set_pair("_reflns_d_resolution_high", f"{min(d):.4f}")
    block.set_pair("_reflns_d_resolution_low", f"{max(d):.4f}")


def _merging_summary(overall):
    """A free text summary of the merging statistics, for a semicolon-delimited
    CIF text field."""

    stats = []

    def add(label, value, fmt):
        if value is not None:
            stats.append(f"{label} {value:{fmt}}")

    add("Rmerge", overall.r_merge, ".4f")
    add("Rmeas", overall.r_meas, ".4f")
    add("CC1/2", overall.cc_one_half, ".4f")
    add("<I/sigma(I)>", overall.i_over_sigma_mean, ".1f")
    add("multiplicity", overall.mean_redundancy, ".1f")
    add("completeness", 100 * overall.completeness, ".1f")
    return "\n".join(
        [
            f"Data reduced with {dials_version()}.",
            f"{overall.n_obs} observations of {overall.n_uniq} unique reflections",
            f"to d_min {overall.d_min:.2f} A.",
            ", ".join(stats) + "%.",
        ]
    )


def _write_absorption(block, experiment_list, params):
    models = [exp.scaling_model for exp in experiment_list]
    if any(m is not None for m in models) and "scale" in params.intensity:
        block.set_pair("_exptl_absorpt_correction_type", "multi-scan")
        components = set()
        for model in models:
            if model is not None:
                components.update(model.components)
        block.set_pair(
            "_exptl_absorpt_process_details",
            gemmi.cif.quote(
                f"{dials_version()} (dials.scale, components: "
                f"{', '.join(sorted(components))})"
            ),
        )
    else:
        block.set_pair("_exptl_absorpt_correction_type", "none")


def _write_reflections(block, reflections, intensities, sigmas, params):
    columns = ["index_h", "index_k", "index_l", "intensity_net", "intensity_u"]
    if params.cif.scale_group_code:
        columns.append("scale_group_code")

    h, k, l = ([str(i) for i in c] for c in zip(*reflections["miller_index"]))
    values = [
        h,
        k,
        l,
        [f"{i:.4f}" for i in intensities],
        [f"{s:.4f}" for s in sigmas],
    ]
    if params.cif.scale_group_code:
        # 1-based, as zero is not a valid SHELX batch number
        values.append([str(i + 1) for i in reflections["id"]])

    loop = block.init_loop("_diffrn_refln_", columns)
    loop.set_all_values(values)

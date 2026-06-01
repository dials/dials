from __future__ import annotations

from iotbx import cif, mtz, pdb

try:
    from mmtbx.command_line.fmodel import (
        fmodel_from_xray_structure_master_params as fmodel_phil,
    )
except ImportError:
    from mmtbx.programs.fmodel import master_phil as fmodel_phil
from mmtbx.utils import fmodel_from_xray_structure

reference_phil_str = """
reference_model {
  k_sol = 0.35
    .type = float
    .help = "Average solvent density to use when calculating the bulk solvent"
            "contribution to the structure factors from a structural model."
            "See Fokine and Urzhumtsev, Acta Cryst. (2002). D58, 1387-1392"
            "for further details on the meaning of this parameter."
    .expert_level = 3
  b_sol = 46.0
    .type = float
    .help = "Average solvent B-factor to use when calculating the bulk solvent"
            "contribution to the structure factors from a structural model."
            "See Fokine and Urzhumtsev, Acta Cryst. (2002). D58, 1387-1392"
            "for further details on the meaning of this parameter."
    .expert_level = 3
}
"""


def intensities_from_reference_file(
    filename, d_min=2.0, wavelength=None, k_sol=0.35, b_sol=46.0
):
    """
    Extract/calculate intensities from a reference file - which may contain data
    or a model, from a pdb, cif or mtz file. For cif, this can be a file
    containing measured data, or a pdb model conforming to mmcif, or an inorgranic
    cif conforming to core_cif.
    """
    if filename.endswith(".pdb"):
        return intensity_array_from_pdb_model_file(
            filename, d_min, wavelength, k_sol, b_sol
        )
    elif filename.endswith(".cif"):
        # need to see if file is a datafile or model file.
        # First try to interpret as a data file (quick and cheap to try)
        try:
            return intensity_array_from_cif_data_file(filename)
        except KeyError:
            return intensity_array_from_cif_model_file(
                filename, d_min, wavelength, k_sol, b_sol
            )
    elif filename.endswith(".mtz"):
        return intensity_array_from_mtz_file(filename)
    else:
        raise ValueError(
            "Unrecognised input format for reference file (expected .pdb, .cif or .mtz)"
        )


def intensities_from_reference_data_file(filename):
    """
    Extract intensities from a reference data file (i.e. cif or mtz containing
    intensities/measured F values.
    """
    if filename.endswith(".cif"):
        return intensity_array_from_cif_data_file(filename)
    elif filename.endswith(".mtz"):
        return intensity_array_from_mtz_file(filename)
    else:
        raise ValueError(
            "Unrecognised input format for reference data file (expected .cif or .mtz)"
        )


def intensities_from_reference_model_file(
    filename, d_min=2.0, wavelength=None, k_sol=0.35, b_sol=46.0
):
    """
    Calculate intensities from a reference model file (i.e. cif or pdb containing
    a crystal structure).
    """
    if filename.endswith(".cif"):
        return intensity_array_from_cif_model_file(
            filename, d_min, wavelength, k_sol, b_sol
        )
    elif filename.endswith(".pdb"):
        return intensity_array_from_pdb_model_file(
            filename, d_min, wavelength, k_sol, b_sol
        )
    else:
        raise ValueError(
            "Unrecognised input format for reference model file (expected .cif or .pdb)"
        )


def intensity_array_from_cif_model_file(
    cif_file, d_min=2.0, wavelength=None, k_sol=0.35, b_sol=46.0
):
    """Return an intensity miller array from a cif file."""
    try:
        # First try to interpret as a MX cif data structure
        xray_structure = pdb.hierarchy.input(cif_file).xray_structure_simple()
    except AssertionError:
        # If this is not understood, try to read as an inorganic cif, using
        # direct methods to calculate the intensities.
        r = cif.reader(file_path=cif_file)
        xray_structure = list(r.build_crystal_structures().values())[0]
        ic = (
            xray_structure.structure_factors(
                anomalous_flag=True, d_min=d_min, algorithm="direct"
            )
            .f_calc()
            .as_intensity_array()
        )
    else:
        # Use the mmtbx methods for quick calculation of fs from model.
        ic = _mmtbx_intensity_from_structure(
            xray_structure, d_min, wavelength, k_sol, b_sol
        )
    return ic


def intensity_array_from_pdb_model_file(
    pdb_file, d_min=2.0, wavelength=None, k_sol=0.35, b_sol=46.0
):
    xray_structure = pdb.hierarchy.input(pdb_file).xray_structure_simple()
    ic = _mmtbx_intensity_from_structure(
        xray_structure, d_min, wavelength, k_sol, b_sol
    )
    return ic


def _mmtbx_intensity_from_structure(
    xray_structure, d_min, wavelength=None, k_sol=0.35, b_sol=46.0
):
    """Use the mmtbx methods for quick calculation of fs from model."""
    if wavelength:
        xray_structure.set_inelastic_form_factors(photon=wavelength, table="sasaki")
    params = fmodel_phil.extract()
    if d_min:
        params.high_resolution = d_min
    params.fmodel.k_sol = k_sol
    params.fmodel.b_sol = b_sol
    if wavelength:
        params.wavelength = wavelength
    fm = fmodel_from_xray_structure(xray_structure, params=params)
    ic = fm.f_model.as_intensity_array()
    return ic


def intensity_array_from_cif_data_file(cif_file):
    f = cif.reader(file_path=cif_file)
    # N.B. in mmcif files, the item is _refln.F_meas (or _refln.F_meas_au),
    # whereas for 'core_cif' inorganic cif the item is _refln_F_meas. This will
    # pick up both cases.
    for ma in f.as_miller_arrays():
        for l in ma.info().labels:
            if "F_meas" in l:
                return ma.as_intensity_array()
    raise KeyError("Unable to find F_meas in cif file")


def intensity_array_from_mtz_file(filename):
    m = mtz.object(filename)
    col_keys = [c.label() for c in m.columns()]
    I_labels_to_try = [
        ["I(+)", "SIGI(+)", "I(-)", "SIGI(-)"],
        ["I", "SIGI"],
        ["IMEAN", "SIGIMEAN"],
        ["IOBS", "SIGIOBS"],  # generated by phenix.cif_as_mtz
    ]
    F_labels_to_try = [
        ["F(+)", "SIGF(+)", "F(-)", "SIGF(-)"],
        ["F", "SIGF"],
        ["FOBS", "SIGFOBS"],  # generated by phenix.cif_as_mtz
    ]
    for labels in I_labels_to_try:
        if all(k in col_keys for k in labels):
            for ma in m.as_miller_arrays():
                if ma.info().labels == labels:
                    return ma
    for labels in F_labels_to_try:
        if all(k in col_keys for k in labels):
            for ma in m.as_miller_arrays():
                if ma.info().labels == labels:
                    return ma.as_intensity_array()
    else:
        raise KeyError("Unable to extract intensity array from mtz file")

from __future__ import annotations

from iotbx import cif, mtz, pdb
from mmtbx.command_line.fmodel import fmodel_from_xray_structure_master_params
from mmtbx.utils import fmodel_from_xray_structure


def intensities_from_reference_file(filename, d_min=2.0, wavelength=None):
    """
    Extract/calculate intensities from a reference file - which may contain data
    or a model, from a pdb, cif or mtz file. For cif, this can be a file
    containing measured data, or a pdb model conforming to mmcif, or an inorgranic
    cif conforming to core_cif.
    """
    if filename.endswith(".pdb"):
        return intensity_array_from_pdb_model_file(filename, d_min, wavelength)
    elif filename.endswith(".cif"):
        # need to see if file is a datafile or model file.
        # First try to interpret as a data file (quick and cheap to try)
        try:
            return intensity_array_from_cif_data_file(filename)
        except KeyError:
            return intensity_array_from_cif_model_file(filename, d_min, wavelength)
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


def intensities_from_reference_model_file(filename, d_min=2.0, wavelength=None):
    """
    Calculate intensities from a reference model file (i.e. cif or pdb containing
    a crystal structure).
    """
    if filename.endswith(".cif"):
        return intensity_array_from_cif_model_file(filename, d_min, wavelength)
    elif filename.endswith(".pdb"):
        return intensity_array_from_pdb_model_file(filename, d_min, wavelength)
    else:
        raise ValueError(
            "Unrecognised input format for reference model file (expected .cif or .pdb)"
        )


def intensity_array_from_cif_model_file(cif_file, d_min=2.0, wavelength=None):
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
        ic = _mmtbx_intensity_from_structure(xray_structure, d_min, wavelength)
    return ic


def intensity_array_from_pdb_model_file(pdb_file, d_min=2.0, wavelength=None):
    xray_structure = pdb.hierarchy.input(pdb_file).xray_structure_simple()
    ic = _mmtbx_intensity_from_structure(xray_structure, d_min, wavelength)
    return ic


def _mmtbx_intensity_from_structure(xray_structure, d_min, wavelength=None):
    """Use the mmtbx methods for quick calculation of fs from model."""
    if wavelength:
        xray_structure.set_inelastic_form_factors(photon=wavelength, table="sasaki")
    params = fmodel_from_xray_structure_master_params.extract()
    params.high_resolution = d_min
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
    cols = m.columns()
    col_dict = {c.label(): c for c in cols}
    # try to make an anomalous intensity array first
    anom_labels = ["I(+)", "SIGI(+)", "I(-)", "SIGI(-)"]
    normal_labels = ["I", "SIGI"]
    mean_labels = ["IMEAN", "SIGIMEAN"]
    F_anom_labels = ["F(+)", "SIGF(+)", "F(-)", "SIGF(-)"]
    F_normal_labels = ["F", "SIGF"]
    I_labels_to_try = [anom_labels, normal_labels, mean_labels]
    F_labels_to_try = [F_anom_labels, F_normal_labels]
    for l in I_labels_to_try:
        if all(k in col_dict.keys() for k in l):
            for ma in m.as_miller_arrays():
                if ma.info().labels == anom_labels:
                    return ma
    for l in F_labels_to_try:
        if all(k in col_dict.keys() for k in l):
            for ma in m.as_miller_arrays():
                if ma.info().labels == anom_labels:
                    return ma.as_intensity_array()
    else:
        raise KeyError("Unable to extract intensity array from mtz file")

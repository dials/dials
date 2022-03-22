"""
Correct integrated intensities to account for attenuation by a diamond anvil cell.

High pressure X-ray diffraction experiments often involve a diamond anvil pressure
cell, in which the sample is sandwiched between two anvils, effectively parallel flat
plates of diamond.  The passage of the incident and diffracted beam through the
anvils results in attenuation of both beams by the diamond by an amount that is
dependent on the path length of each beam through each anvil.

This utility calculates these path lengths and boosts the integrated reflection
intensities to remove the calculated effect of the diamond attenuation.

It is intended that this program be used to correct reflection intensities after
integration but before scaling.  Call it on the output of dials.integrate.

Examples::

  dials.anvil_correction integrated.expt integrated.refl

  dials.anvil_correction integrated.expt integrated.refl thickness=1.2 normal=1,0,0
"""


from __future__ import annotations

import logging
import sys
from typing import List, Sequence, SupportsFloat

import numpy as np
from scipy.spatial.transform import Rotation

import libtbx.phil
from cctbx.eltbx import attenuation_coefficient
from dxtbx.model import Experiment

import dials.util.log
from dials.array_family import flex
from dials.util.multi_dataset_handling import (
    parse_multiple_datasets,
    sort_tables_to_experiments_order,
)
from dials.util.options import ArgumentParser, flatten_experiments, flatten_reflections

Vector = Sequence[SupportsFloat]


logger = logging.getLogger("dials.command_line.anvil_correction")

phil_scope = libtbx.phil.parse(
    """
    anvil
        .caption = "Properties of the mounted diamond anvils"
    {
        density = 3510
            .type = float
            .help = "The density of the anvil material in kg per cubic metre.  " \
                    "The default is the typical density of synthetic diamond."

        thickness = 1.5925
            .type = float
            .help = "The thickness in mm of each anvil in the pressure cell.  " \
                    "The default is the thickness of the pressure cells in use on " \
                    "beamline I19 at Diamond Light Source."

        normal = 0, 1, 0
            .type = floats(size=3)
            .help = "A 3-vector orthogonal to the anvil surfaces in the laboratory " \
                    "frame when the goniometer is at zero datum, i.e. the axes are " \
                    "all at zero degrees.  The vector may be given un-normalised."
    }

    output {
        experiments = None
            .type = path
            .help = "The output experiment list file name. If None, don't output an " \
                    "experiment list file."
        reflections = corrected.refl
            .type = path
            .help = "The output reflection table file."
        log = dials.anvil_correction.log
            .type = path
    }
    """
)

# Get the tabulated NIST mass attenuation coefficient data for carbon.
carbon_attenuation_data = attenuation_coefficient.get_table("C")


def goniometer_rotation(
    experiment: Experiment, reflections: flex.reflection_table
) -> Rotation:
    """
    Calculate the goniometer rotation operator for each reflection.

    For each reflection, find the rotation operator that describes the position in
    the lab frame that the sample was in when the reflection was measured.

    Following the DXTBX model of a goniometer, whereby a scan is only possible
    around one physical axis at a time, the goniostat rotation operator can be
    represented as
      R = S ∘ R' ∘ F.
    Here:
        * S is the static 'setting rotation', the operator denoting the position of all
        parent axes of the scan axis, which hence defines the orientation of the scan
        axis;
        * R' is the the operator denoting the scan rotation.  It has a different
        value for each reflection, recording the scan position corresponding to each
        reflection centroid;
        * F is the static 'fixed rotation', denoting the orientation of all child axes
        of the scan axis.

    Args:
        experiment:  The DXTBX experiment object corresponding to the scan.
        reflections:  A table of reflections at which to calculate the rotations.

    Returns:
        An array of rotation operators, one per reflection in the reflection table.
    """
    # Get the axis of the scan rotation.
    rotation_axis = experiment.goniometer.get_rotation_axis_datum()
    # For each reflection, get the angle of the scan rotation.
    angles = reflections["xyzobs.mm.value"].parts()[2]
    # Construct a rotation vector (parallel with the rotation axis and with
    # magnitude equal to the rotation angle) for each reflection.
    # The shape of this array is (N, 3), where N is the number of reflections.
    rotvecs = np.outer(angles, rotation_axis)
    # Create a rotation operator for each scan rotation (i.e. one per reflection).
    # In the notation of dxtbx/model/goniometer.h, this is R.
    scan_rotation = Rotation.from_rotvec(rotvecs)

    # Get the setting rotation.
    # In the notation of dxtbx/model/goniometer.h, this is S.
    set_rotation = np.array(experiment.goniometer.get_setting_rotation()).reshape(3, 3)
    set_rotation = Rotation.from_matrix(set_rotation)

    # Create a rotation operator for those axes that are fixed throughout the scan.
    # In the notation of dxtbx/model/goniometer.h, this is F.
    fixed_rotation = np.array(experiment.goniometer.get_fixed_rotation()).reshape(3, 3)
    fixed_rotation = Rotation.from_matrix(fixed_rotation)

    # Calculate the rotation operator representing the goniometer orientation for each
    # reflection.  In the notation of dxtbx/model/goniometer.h this is S ∘ R ∘ F.
    return set_rotation * scan_rotation * fixed_rotation


def attenuation_correction(
    experiment: Experiment,
    reflections: flex.reflection_table,
    dac_norm: Vector,
    thickness: float,
    density: float,
) -> flex.double:
    """
    Calculate the correction factors for attenuation by a diamond anvil cell.

    Take an experiment object and reflection table containing integrated but unscaled
    diffraction data, estimated the factor by which each of the integrated
    reflection intensities must be boosted to correct for attenuation of the incident
    and diffracted beams in the diamond anvils.

    Args:
        experiment:  An experiment from integrated data.
        reflections:  A table of integrated reflections.
        dac_norm:  A 3-vector representing the normal to the anvil surfaces in the
                   laboratory frame when the goniometer is at zero datum, i.e. the axes
                   are all at 0°.  The vector is assumed to be normalised.
        thickness:  The thickness of each diamond anvil (assumed equal).
        density:  The density of the anvil material in g·cm⁻³
                  (units chosen for consistency with the NIST tables of X-ray mass
                  attenuation coefficients, https://dx.doi.org/10.18434/T4D01F).

    Returns:
        An array of correction factors for the integrated reflection intensities.
    """
    # Get the wavelength.
    wavelength = experiment.beam.get_wavelength()
    # Get the mass attenuation coefficient in cm²·g⁻¹ from the NIST tables.
    mass_atten_coeff = carbon_attenuation_data.mu_rho_at_angstrom(wavelength)
    # Get the linear attenuation coefficient in mm⁻¹.
    linear_atten_coeff = density * mass_atten_coeff / 10  # mm⁻¹

    # Find the orientation of the diamond anvil cell for each reflection.  Since the
    # cell is fixed to the goniometer, this is just S × R × F × ̠̂n.
    # This is an array of shape (N, 3), i.e. a 3-vector for each reflection.
    dac_orientation = goniometer_rotation(experiment, reflections).apply(dac_norm)

    # Get the normalised incident and diffracted beam vectors ̠̂s₀ and ̠̂s₁.
    # Naturally, there will be one vale of ̠̂s₁ for each reflection.
    s0_norm = experiment.beam.get_unit_s0()  # An array of shape (3).
    s1_norm = reflections["s1"] * wavelength  # Shape (N, 3).

    # Get the scalar product of the diamond anvil cell orientation with each of
    # ̠̂s₀ & ̠̂s₁ for each reflection.
    incident_cosine = np.dot(dac_orientation, s0_norm)  # An array of shape (N).
    # Dot product of each of an array of N 3-vectors with its counterpart from another
    # array of N 3-vectors.  Shape (N, 3) • (N, 3) → N.
    # Equivalent to np.sum(dac_orientation * s1_norm, axis=1), but quicker.
    diffracted_cosine = np.einsum("ij,ij->i", dac_orientation, s1_norm)  # Shape (N)

    # Get the path length through the anvil of the incident and reflected beams.
    l0 = thickness / np.abs(incident_cosine)  # mm.  An array of shape (N).
    l1 = thickness / np.abs(diffracted_cosine)  # mm.  An array of shape (N).

    # Derive the factor by which we estimate the intensities to have been attenuated
    # by the anvils.  This is an array of shape (N).
    l_tot = l0 + l1
    return flex.double(np.exp(linear_atten_coeff * l_tot))


def correct_intensities_for_dac_attenuation(
    experiment: Experiment,
    reflections: flex.reflection_table,
    dac_norm: Vector,
    thickness: float,
    density: float = 3.51,
) -> None:
    """
    Boost integrated intensities to account for attenuation by a diamond anvil cell.

    Take an experiment object and reflection table containing integrated but unscaled
    diffraction data, boost the integrated reflection intensities to correct for the
    estimated attenuation due to the passage of the incident and diffracted beams
    through the diamond anvils.

    Args:
        experiment:  An experiment from integrated data.
        reflections:  A table of reflections from integrated data.
        dac_norm:  A 3-vector representing the normal to the anvil surfaces in the
                   laboratory frame when the goniometer is at zero datum, i.e. the axes
                   are all at 0°.  The vector is assumed to be normalised.
        thickness:  The thickness of each diamond anvil (assumed equal).
        density:  The density of the anvil material in g·cm⁻³
                  (units chosen for consistency with the NIST tables of X-ray mass
                  attenuation coefficients, https://dx.doi.org/10.18434/T4D01F).
                  Defaults to a typical value for synthetic diamond.
    """
    # Select only those reflections whose intensities after integration ought to be
    # meaningful.
    sel = reflections.get_flags(reflections.flags.integrated, all=False)
    sel = sel.iselection()
    refls_sel = reflections.select(sel)

    correction = attenuation_correction(
        experiment, refls_sel, dac_norm, thickness, density
    )

    # We need only correct non-null values for each integration method.
    prf_subsel = refls_sel.get_flags(refls_sel.flags.integrated_prf)
    sum_subsel = refls_sel.get_flags(refls_sel.flags.integrated_sum)

    # Correct the measured intensities and variances for this attenuation.
    methods = {"prf": prf_subsel, "sum": sum_subsel}
    quantities = {"value": correction, "variance": flex.pow2(correction)}
    for method, subsel in methods.items():
        setting_subsel = sel.select(subsel)
        for quantity, factor in quantities.items():
            col = f"intensity.{method}.{quantity}"
            corrected = (refls_sel[col] * factor).select(subsel)
            try:
                reflections[col].set_selected(setting_subsel, corrected)
            except KeyError:
                pass


@dials.util.show_mail_handle_errors()
def run(args: List[str] = None, phil: libtbx.phil.scope = phil_scope) -> None:
    """
    Run dials.anvil_correction as from the command line.

    Take integrated experiment lists and reflection tables and correct the all the
    integrated intensities for the estimated attenuation by the diamond anvils.

    Args:
        args: The arguments supplied by the user (default: sys.argv[1:]).
        phil: The PHIL scope definition (default: phil_scope, the master PHIL scope
              for this program).
    """
    usage = "dials.anvil_correction [options] integrated.expt integrated.refl"

    parser = ArgumentParser(
        usage=usage,
        phil=phil,
        read_reflections=True,
        read_experiments=True,
        check_format=False,
        epilog=__doc__,
    )

    params, options = parser.parse_args(args=args, show_diff_phil=False)

    # Log the difference between the PHIL scope definition and the active PHIL scope,
    # which will include the parsed user inputs.
    diff_phil = parser.diff_phil.as_str()
    if diff_phil:
        logger.info("The following parameters have been modified:\n%s", diff_phil)

    # Check that at least one reflection table and experiment list have been provided.
    input_errors = []
    if not params.input.experiments:
        input_errors.append(
            "Please provide at least one valid experiment list (.expt) file."
        )
    if not params.input.reflections:
        input_errors.append(
            "Please provide at least one valid reflection table (.refl) file."
        )
    if input_errors:
        sys.exit("\n".join([parser.format_help()] + input_errors))

    if not np.linalg.norm(params.anvil.normal):
        sys.exit("It seems you have provided a surface normal vector with zero length.")

    # Check that the anvil surface normal really is normalised.
    dac_norm = params.anvil.normal / np.linalg.norm(params.anvil.normal)

    # Configure the logging.
    dials.util.log.config(options.verbose, logfile=params.output.log)

    # These functions are commonly used to collate the input.
    experiments = flatten_experiments(params.input.experiments)
    reflections_list = flatten_reflections(params.input.reflections)
    # Work around parse_multiple_datasets dropping unindexed reflections.
    unindexed = flex.reflection_table()
    for r_table in reflections_list:
        unindexed.extend(r_table.select(r_table["id"] == -1))
    # Get a single reflection table per experiment object.
    reflections_list = parse_multiple_datasets(reflections_list)
    reflections_list = sort_tables_to_experiments_order(reflections_list, experiments)

    # Record the density of diamond in g·cm⁻³ (for consistency with NIST tables,
    # https://doi.org/10.18434/T4D01F).
    density = params.anvil.density / 1000  # g·cm⁻³

    # Correct for the attenuation of the incident and diffracted beams by the diamond
    # anvil cell.
    logger.info(
        "Correcting integrated reflection intensities for attenuation by the diamond "
        "anvil cell."
    )
    for experiment, reflections in zip(experiments, reflections_list):
        correct_intensities_for_dac_attenuation(
            experiment, reflections, dac_norm, params.anvil.thickness, density
        )
    logger.info("Done.")

    # Do optional experiment list file output here.
    if params.output.experiments:
        logger.info("Writing the experiment list to %s", params.output.experiments)
        experiments.as_file(params.output.experiments)
    logger.info("Writing the reflection table to %s", params.output.reflections)
    # Collate reflections into a single reflection table and save it to file.
    reflections = unindexed
    for r_table in reflections_list:
        reflections.extend(r_table)
    del reflections_list
    reflections.as_file(params.output.reflections)


# Keep this minimal.  Try to keep the command-line behaviour neatly encapsulated in run.
if __name__ == "__main__":
    run()

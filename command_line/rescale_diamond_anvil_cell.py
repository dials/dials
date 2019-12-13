# coding: utf-8

"""
Boost integrated intensities to account for attenuation by a diamond anvil cell.

High pressure X-ray diffraction experiments often involve a diamond anvil pressure
cell, in which the sample is sandwiched between two anvils, effectively parallel flat
plates of diamond.  The passage of the incident and diffracted beam through the
anvils results in attenuation of both beams by the diamond by an amount that is
dependent on the path length of each beam through each anvil.

This utility calculates these path lengths and boosts the integrated reflection
intensities to remove the calculated effect of the diamond attenuation.

It is intended that this program be used to correct reflection intensities after
integration but before scaling.  Call it on the output of dials.integrate.
"""

from __future__ import absolute_import, division, print_function

import logging
import sys

import numpy as np
from scipy.spatial.transform import Rotation

from cctbx.eltbx import attenuation_coefficient
from dials.array_family import flex
import dials.util.log
from dials.util.options import flatten_experiments, flatten_reflections, OptionParser
from dxtbx.model import ExperimentList
import libtbx.phil

try:
    from typing import List, Sequence, SupportsFloat
except ImportError:
    pass
else:
    Vector = Sequence[SupportsFloat]


logger = logging.getLogger("dials.rescale_diamond_anvil_cell")

phil_scope = libtbx.phil.parse(
    u"""
    anvil
        .caption = "Properties of the mounted diamond anvils"
    {
        density = 3510
            .type = float
            .help = "The density of the anvil material in kg per cubic metre.  "
                    "The default is the typical density of diamond."

        thickness = 1.5925
            .type = float
            .help = "The thickness in mm of each anvil in the pressure cell.\n"
                    "The default is the thickness of the pressure cells in use on "
                    "beamline I19 at Diamond Light Source."

        normal = 0, 1, 0
            .type = floats(size=3)
            .help = "A 3-vector representing the normal to the anvil surfaces in the "
                    "laboratory frame when the goniometer is at zero datum, i.e. the "
                    "axes are all at zero degrees.  The vector may be given "
                    "un-normalised."
    }

    output {
        experiments = corrected.expt
            .type = path
            .help = "The output experiment list file name.\n"
                    "If None, don't output an experiment list file."
        reflections = corrected.refl
            .type = path
            .help = "The output reflection table file."
        log = dials.command_name.log
            .type = path
    }
    """
)

# Get the tabulated NIST mass attenuation coefficient data for carbon.
carbon_attenuation_data = attenuation_coefficient.get_table("C")


def correct_intensities_for_dac_attenuation(
    elist,  # type: ExperimentList
    rtable,  # type: flex.reflection_table
    dac_norm,  # type: Vector
    thickness,  # type: float
    density=3.51,  # type: float  # g·cm⁻³
):  # type: (...) -> None
    u"""
    Boost integrated intensities to account for attenuation by a diamond anvil cell.

    Take an experiment list and reflection table containing integrated but unscaled
    diffraction data, boost the integrated reflection intensities to correct for the
    estimated attenuation due to the passage of the incident and diffracted beams
    through the diamond anvils.

    Args:
        elist:  An experiment list from integrated data.
        rtable:  A table of integrated reflections.
        dac_norm:  A 3-vector representing the normal to the anvil surfaces in the
                   laboratory frame when the goniometer is at zero datum, i.e. the axes
                   are all at 0°.  The vector is assumed to be normalised.
        thickness:  The thickness of each diamond anvil (assumed equal).
        density:  The density of the anvil material in g·cm⁻³
                  (units chosen for consistency with the NIST tables of X-ray mass
                  attenuation coefficients, https://dx.doi.org/10.18434/T4D01F).
                  Defaults to a typical value for synthetic diamond.
    """
    # Iterate over all experiments in the experiment list.
    for i, expt in enumerate(elist):
        # Get the wavelength.
        wavelength = expt.beam.get_wavelength()
        # Get the mass attenuation coefficient in cm²·g⁻¹ from the NIST tables.
        mass_atten_coeff = carbon_attenuation_data.mu_rho_at_angstrom(wavelength)
        # Get the linear attenuation coefficient in mm⁻¹.
        linear_atten_coeff = density * mass_atten_coeff / 10  # mm⁻¹

        # Select the reflections with the correct experiment ID and only those whose
        # intensities after integration ought to be meaningful.
        sel = rtable["id"] == i
        sel &= rtable.get_flags(rtable.flags.integrated, all=False)
        sel = sel.iselection()
        refls = rtable.select(sel)

        # Get the setting rotation.
        # In the notation of dxtbx/model/goniometer.h, this is S.
        set_rotation = np.array(expt.goniometer.get_setting_rotation()).reshape(3, 3)
        set_rotation = Rotation.from_dcm(set_rotation)

        # Get the axis of the scan rotation.
        rotation_axis = expt.goniometer.get_rotation_axis_datum()
        # For each reflection, get the angle of the scan rotation.
        angles = refls["xyzobs.mm.value"].parts()[2]
        # Construct a rotation vector (parallel with the rotation axis and with
        # magnitude equal to the rotation angle) for each reflection.
        # The shape of this array is (N, 3), where N is the number of reflections.
        rotvecs = np.outer(angles, rotation_axis)
        # Remove redundant things that scale with N.
        del angles
        # Create a rotation operator for each scan rotation (i.e. one per reflection).
        # In the notation of dxtbx/model/goniometer.h, this is R.
        scan_rotation = Rotation.from_rotvec(rotvecs)
        # Remove redundant things that scale with N.
        del rotvecs

        # Create a rotation operator for those axes that are fixed throughout the scan.
        # In the notation of dxtbx/model/goniometer.h, this is F.
        fixed_rotation = np.array(expt.goniometer.get_fixed_rotation()).reshape(3, 3)
        fixed_rotation = Rotation.from_dcm(fixed_rotation)

        # Get the rotation operator representing the goniometer orientation for each
        # reflection.  In the notation of dxtbx/model/goniometer.h this is S × R × F.
        rotation = set_rotation * scan_rotation * fixed_rotation
        # Remove redundant things that scale with N.
        del scan_rotation

        # Find the orientation of the diamond anvil cell for each reflection.  Since the
        # cell is fixed to the goniometer, this is just S × R × F × ̠̂n.
        # This is an array of shape (N, 3), i.e. a 3-vector for each reflection.
        dac_orientation = rotation.apply(dac_norm)
        # Remove redundant things that scale with N.
        del rotation

        # Get the normalised incident and diffracted beam vectors ̠̂s₀ and ̠̂s₁.
        # Naturally, there will be one vale of ̠̂s₁ for each reflection.
        s0_norm = expt.beam.get_unit_s0()  # An array of shape (3).
        s1_norm = refls["s1"] * wavelength  # Shape (N, 3).

        # Get the scalar product of the diamond anvil cell orientation with each of
        # ̠̂s₀ & ̠̂s₁ for each reflection.
        incident_cosine = np.dot(dac_orientation, s0_norm)  # An array of shape (N).
        diffracted_cosine = np.einsum("ij,ij->i", dac_orientation, s1_norm)  # Shape (N)
        # Remove redundant things that scale with N.
        del s1_norm, dac_orientation

        # Get the path length through the anvil of the incident and reflected beams.
        l0 = thickness / np.abs(incident_cosine)  # mm.  An array of shape (N).
        # Remove redundant things that scale with N.
        del incident_cosine
        l1 = thickness / np.abs(diffracted_cosine)  # mm.  An array of shape (N).
        # Remove redundant things that scale with N.
        del diffracted_cosine

        # Derive the factor by which we estimate the intensities to have been attenuated
        # by the anvils.  This is an array of shape (N).
        transmission = np.exp(-linear_atten_coeff * l0)
        transmission *= np.exp(-linear_atten_coeff * l1)
        # Remove redundant things that scale with N.
        del l0, l1

        # Correct the measured intensities for this attenuation.
        for col in "intensity.prf.value", "intensity.sum.value":
            try:
                rtable[col].set_selected(sel, refls[col] / flex.double(transmission))
            except KeyError:
                pass
        # Correct the measured variances accordingly.
        for col in "intensity.prf.variance", "intensity.sum.variance":
            try:
                rtable[col].set_selected(
                    sel, refls[col] / flex.double(np.square(transmission))
                )
            except KeyError:
                pass


def run(args=None, phil=phil_scope):  # type: (List[str], libtbx.phil.scope) -> None
    """
    Run dev.dials.rescale_diamond_anvil_cell as from the command line.

    Take integrated experiment lists and reflection tables and correct the all the
    integrated intensities for the estimated attenuation by the diamond anvils.

    Args:
        args: The arguments supplied by the user (default: sys.argv[1:]).
        phil: The PHIL scope definition (default: phil_scope, the master PHIL scope
              for this program).
    """
    usage = (
        "dials.rescale_diamond_anvil_cell [options] integrated.expt " "integrated.refl"
    )

    parser = OptionParser(
        usage=usage,
        phil=phil,
        read_reflections=True,
        read_experiments=True,
        check_format=False,
        epilog=__doc__,
    )

    params, options = parser.parse_args(args=args, show_diff_phil=False)

    # Configure the logging.
    dials.util.log.config(options.verbose, logfile=params.output.log)

    # Log the difference between the PHIL scope definition and the active PHIL scope,
    # which will include the parsed user inputs.
    diff_phil = parser.diff_phil.as_str()
    if diff_phil:
        logger.info("The following parameters have been modified:\n%s", diff_phil)

    # These functions are commonly used to collate the input.
    experiments = flatten_experiments(params.input.experiments)
    reflections_list = flatten_reflections(params.input.reflections)
    reflections = reflections_list[0]
    for r_table in reflections_list[1:]:
        reflections.extend(r_table)
    del reflections_list

    # Check that the anvil surface normal really is normalised.
    try:
        dac_norm = params.anvil.normal / np.linalg.norm(params.anvil.normal)
    except ZeroDivisionError:
        sys.exit("It seems you have provided a surface normal vector with zero length.")

    # Record the density of diamond in g·cm⁻³ (for consistency with NIST tables,
    # https://dx.doi.org/10.18434/T4D01F).
    density = params.anvil.density / 1000  # g·cm⁻³

    # Correct for the attenuation of the incident and diffracted beams by the diamond
    # anvil cell.
    logger.info(
        "Correcting integrated reflection intensities for attenuation by the diamond "
        "anvil cell."
    )
    correct_intensities_for_dac_attenuation(
        experiments, reflections, dac_norm, params.anvil.thickness, density
    )
    logger.info("Done.")

    # Do the file output here.
    if params.output.experiments:
        logger.info("Writing the experiment list to %s", params.output.experiments)
        experiments.as_file(params.output.experiments)
    logger.info("Writing the reflection table to %s", params.output.reflections)
    reflections.as_file(params.output.reflections)


# Keep this minimal.  Try to keep the command-line behaviour neatly encapsulated in run.
if __name__ == "__main__":
    with dials.util.show_mail_on_error():
        run()

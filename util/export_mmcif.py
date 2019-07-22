from __future__ import absolute_import, division, print_function

import datetime
import logging
import math

import dials.util.version
from dials.util.filter_reflections import filter_reflection_table
import iotbx.cif.model
from cctbx.sgtbx import bravais_types
from cctbx import miller
from cctbx import crystal as cctbxcrystal
from iotbx.merging_statistics import dataset_statistics
from libtbx import Auto

logger = logging.getLogger(__name__)
RAD2DEG = 180.0 / math.pi


class MMCIFOutputFile(object):
    """
    Class to output experiments and reflections as MMCIF file
    """

    def __init__(self, params):
        """
        Init with the filename
        """
        self._cif = iotbx.cif.model.cif()
        self.params = params

    def write(self, experiments, reflections):
        """
        Write the experiments and reflections to file
        """

        # if mmmcif filename is auto, then choose scaled.cif or integrated.cif
        if self.params.mmcif.hklout in (None, Auto, "auto"):
            if ("intensity.scale.value" in reflections) and (
                "intensity.scale.variance" in reflections
            ):
                filename = "scaled.cif"
                logger.info(
                    "Data appears to be scaled, setting mmcif.hklout = 'scaled_unmerged.cif'"
                )
            else:
                filename = "integrated.cif"
                logger.info(
                    "Data appears to be unscaled, setting mmcif.hklout = 'integrated.cif'"
                )

        # Select reflections
        selection = reflections.get_flags(reflections.flags.integrated, all=True)
        reflections = reflections.select(selection)

        # Filter out bad variances and other issues, but don't filter on ice rings
        # or alter partialities.

        ### Assumes you want to apply the lp and dqe corrections to sum and prf
        ### Do we want to combine partials?
        reflections = filter_reflection_table(
            reflections,
            self.params.intensity,
            combine_partials=False,
            partiality_threshold=0.0,
            d_min=self.params.mtz.d_min,
        )

        # Get the cif block
        cif_block = iotbx.cif.model.block()

        # Audit trail
        dials_version = dials.util.version.dials_version()
        cif_block["_audit.creation_method"] = dials_version
        cif_block["_audit.creation_date"] = datetime.date.today().isoformat()
        cif_block["_computing.data_reduction"] = (
            "%s (Winter, G. et al., 2018)" % dials_version
        )
        cif_block[
            "_publ.section_references"
        ] = "Winter, G. et al. (2018) Acta Cryst. D74, 85-97."

        # Hard coding X-ray
        cif_block["_pdbx_diffrn_data_section.id"] = "dials"
        cif_block["_pdbx_diffrn_data_section.type_scattering"] = "x-ray"
        cif_block["_pdbx_diffrn_data_section.type_merged"] = "false"
        cif_block["_pdbx_diffrn_data_section.type_scaled"] = str(
            "scale" in self.params.intensity
        ).lower()

        # FIXME Haven't put in any of these bits yet
        #
        #  Facility/beamline proposal tracking details
        #
        # cif_block["_pdbx_diffrn_data_section_experiment.ordinal"] = 1
        # cif_block["_pdbx_diffrn_data_section_experiment.data_section_id"] = "dials"
        # cif_block["_pdbx_diffrn_data_section_experiment.proposal_id"] = "<PROPOSAL ID>

        # Facility/beamline details for this data collection
        #
        # cif_block["_pdbx_diffrn_data_section_site.data_section_id"] = 'dials'
        # cif_block["_pdbx_diffrn_data_section_site.facility"] = "DIAMOND"
        # cif_block["_pdbx_diffrn_data_section_site.beamline"] = "VMX-M"
        # cif_block["_pdbx_diffrn_data_section_site.collection_date"] = scan.epochs()[0]
        # cif_block["_pdbx_diffrn_data_section_site.detector"] = detector[0].name()
        # cif_block["_pdbx_diffrn_data_section_site.detector_type"] = detector[0].type()

        # Write the crystal information
        cif_loop = iotbx.cif.model.loop(
            header=(
                "_pdbx_diffrn_unmerged_cell.ordinal",
                "_pdbx_diffrn_unmerged_cell.crystal_id",
                "_pdbx_diffrn_unmerged_cell.wavelength",
                "_pdbx_diffrn_unmerged_cell.cell_length_a",
                "_pdbx_diffrn_unmerged_cell.cell_length_b",
                "_pdbx_diffrn_unmerged_cell.cell_length_c",
                "_pdbx_diffrn_unmerged_cell.cell_angle_alpha",
                "_pdbx_diffrn_unmerged_cell.cell_angle_beta",
                "_pdbx_diffrn_unmerged_cell.cell_angle_gamma",
                "_pdbx_diffrn_unmerged_cell.Bravais_lattice",
            )
        )
        crystals = experiments.crystals()
        crystal_to_id = {crystal: i + 1 for i, crystal in enumerate(crystals)}
        for i, exp in enumerate(experiments):
            crystal = exp.crystal
            crystal_id = crystal_to_id[crystal]
            wavelength = exp.beam.get_wavelength()
            a, b, c, alpha, beta, gamma = crystal.get_unit_cell().parameters()
            latt_type = str(
                bravais_types.bravais_lattice(group=crystal.get_space_group())
            )
            cif_loop.add_row(
                (i + 1, crystal_id, wavelength, a, b, c, alpha, beta, gamma, latt_type)
            )
            cif_block.add_loop(cif_loop)

        # Write the scan information
        cif_loop = iotbx.cif.model.loop(
            header=(
                "_pdbx_diffrn_scan.scan_id",
                "_pdbx_diffrn_scan.crystal_id",
                "_pdbx_diffrn_scan.image_id_begin",
                "_pdbx_diffrn_scan.image_id_end",
                "_pdbx_diffrn_scan.scan_angle_begin",
                "_pdbx_diffrn_scan.scan_angle_end",
            )
        )
        for i, exp in enumerate(experiments):
            scan = exp.scan
            crystal_id = crystal_to_id[exp.crystal]
            image_range = scan.get_image_range()
            osc_range = scan.get_oscillation_range(deg=True)
            cif_loop.add_row(
                (
                    i + 1,
                    crystal_id,
                    image_range[0],
                    image_range[1],
                    osc_range[0],
                    osc_range[1],
                )
            )
            cif_block.add_loop(cif_loop)

        # Make a dict of unit_cell parameters
        unit_cell_parameters = {}
        if crystal.num_scan_points > 1:
            for i in range(crystal.num_scan_points):
                a, b, c, alpha, beta, gamma = crystal.get_unit_cell_at_scan_point(
                    i
                ).parameters()
                unit_cell_parameters[i] = (a, b, c, alpha, beta, gamma)
        else:
            unit_cell_parameters[0] = (a, b, c, alpha, beta, gamma)

        ### _pdbx_diffrn_image_proc has been removed from the dictionary extension.
        ### Keeping this section commented out as it may be added back in some
        ### form in future
        #
        # Write the image data
        # scan = experiments[0].scan
        # z0 = scan.get_image_range()[0]
        #
        # cif_loop = iotbx.cif.model.loop(
        #  header=("_pdbx_diffrn_image_proc.image_id",
        #          "_pdbx_diffrn_image_proc.crystal_id",
        #          "_pdbx_diffrn_image_proc.image_number",
        #          "_pdbx_diffrn_image_proc.phi_value",
        #          "_pdbx_diffrn_image_proc.wavelength",
        #          "_pdbx_diffrn_image_proc.cell_length_a",
        #          "_pdbx_diffrn_image_proc.cell_length_b",
        #          "_pdbx_diffrn_image_proc.cell_length_c",
        #          "_pdbx_diffrn_image_proc.cell_angle_alpha",
        #          "_pdbx_diffrn_image_proc.cell_angle_beta",
        #          "_pdbx_diffrn_image_proc.cell_angle_gamma"))
        # for i in range(len(scan)):
        #  z = z0 + i
        #  if crystal.num_scan_points > 1:
        #    a, b, c, alpha, beta, gamma = unit_cell_parameters[i]
        #  else:
        #    a, b, c, alpha, beta, gamma = unit_cell_parameters[0]
        #  # phi is the angle at the image centre
        #  phi = scan.get_angle_from_image_index(z + 0.5, deg=True)
        #  cif_loop.add_row((i+1, 1, z, phi, wavelength,
        #                    a, b, c, alpha, beta, gamma))
        # cif_block.add_loop(cif_loop)

        # Write reflection data
        # Required columns
        header = (
            "_pdbx_diffrn_unmerged_refln.reflection_id",
            "_pdbx_diffrn_unmerged_refln.scan_id",
            "_pdbx_diffrn_unmerged_refln.image_id_begin",
            "_pdbx_diffrn_unmerged_refln.image_id_end",
            "_pdbx_diffrn_unmerged_refln.index_h",
            "_pdbx_diffrn_unmerged_refln.index_k",
            "_pdbx_diffrn_unmerged_refln.index_l",
        )

        headernames = {
            "scales": "_pdbx_diffrn_unmerged_refln.scale_value",
            "intensity.scale.value": "_pdbx_diffrn_unmerged_refln.intensity_meas",
            "intensity.scale.sigma": "_pdbx_diffrn_unmerged_refln.intensity_sigma",
            "intensity.sum.value": "_pdbx_diffrn_unmerged_refln.intensity_sum",
            "intensity.sum.sigma": "_pdbx_diffrn_unmerged_refln.intensity_sum_sigma",
            "intensity.prf.value": "_pdbx_diffrn_unmerged_refln.intensity_prf",
            "intensity.prf.sigma": "_pdbx_diffrn_unmerged_refln.intensity_prf_sigma",
            "angle": "_pdbx_diffrn_unmerged_refln.scan_angle_reflection",
            "partiality": "_pdbx_diffrn_unmerged_refln.partiality",
        }

        variables_present = []
        if "scale" in self.params.intensity:
            reflections["scales"] = 1.0 / reflections["inverse_scale_factor"]
            reflections["intensity.scale.sigma"] = (
                reflections["intensity.scale.variance"] ** 0.5
            )
            variables_present.extend(
                ["scales", "intensity.scale.value", "intensity.scale.sigma"]
            )
        if "sum" in self.params.intensity:
            reflections["intensity.sum.sigma"] = (
                reflections["intensity.sum.variance"] ** 0.5
            )
            variables_present.extend(["intensity.sum.value", "intensity.sum.sigma"])
        if "profile" in self.params.intensity:
            reflections["intensity.prf.sigma"] = (
                reflections["intensity.prf.variance"] ** 0.5
            )
            variables_present.extend(["intensity.prf.value", "intensity.prf.sigma"])

        # Should always exist
        reflections["angle"] = reflections["xyzcal.mm"].parts()[2] * RAD2DEG
        variables_present.extend(["angle"])

        if "partiality" in reflections:
            variables_present.extend(["partiality"])

        for name in variables_present:
            if name in reflections:
                header += (headernames[name],)

        if "scale" in self.params.intensity:
            # Write dataset_statistics - first make a miller array
            crystal_symmetry = cctbxcrystal.symmetry(
                space_group=experiments[0].crystal.get_space_group(),
                unit_cell=experiments[0].crystal.get_unit_cell(),
            )
            miller_set = miller.set(
                crystal_symmetry=crystal_symmetry,
                indices=reflections["miller_index"],
                anomalous_flag=False,
            )
            i_obs = miller.array(miller_set, data=reflections["intensity.scale.value"])
            i_obs.set_observation_type_xray_intensity()
            i_obs.set_sigmas(reflections["intensity.scale.sigma"])
            i_obs.set_info(
                miller.array_info(source="DIALS", source_type="reflection_tables")
            )

            result = dataset_statistics(
                i_obs=i_obs,
                crystal_symmetry=crystal_symmetry,
                use_internal_variance=False,
                eliminate_sys_absent=False,
            )

            cif_block.update(result.as_cif_block())

        cif_loop = iotbx.cif.model.loop(header=header)

        for i, r in enumerate(reflections):
            refl_id = i + 1
            scan_id = r["id"] + 1
            _, _, _, _, z0, z1 = r["bbox"]
            h, k, l = r["miller_index"]
            variable_values = tuple((r[name]) for name in variables_present)
            cif_loop.add_row((refl_id, scan_id, z0, z1, h, k, l) + variable_values)
        cif_block.add_loop(cif_loop)

        # Add the block
        self._cif["dials"] = cif_block

        # Print to file
        with open(filename, "w") as fh:
            self._cif.show(out=fh)

        # Log
        logger.info("Wrote reflections to %s" % filename)

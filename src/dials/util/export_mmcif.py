from __future__ import annotations

import bz2
import datetime
import gzip
import logging
import lzma
import math
import time

import iotbx.cif.model
from cctbx import crystal as cctbxcrystal
from cctbx import miller
from cctbx.sgtbx import bravais_types
from iotbx.merging_statistics import dataset_statistics
from libtbx import Auto
from scitbx.array_family import flex

import dials.util.version
from dials.util.filter_reflections import filter_reflection_table

logger = logging.getLogger(__name__)
RAD2DEG = 180.0 / math.pi


class MMCIFOutputFile:
    """
    Class to output experiments and reflections as MMCIF file
    """

    def __init__(self, params):
        """
        Init with the parameters
        """
        self._cif = iotbx.cif.model.cif()
        self.params = params
        self._fmt = "%6i %2i %5i %5i %-2i %-2i %-2i"

    def write(self, experiments, reflections):
        """
        Write the experiments and reflections to file
        """

        # if mmcif filename is auto, then choose scaled.cif or integrated.cif
        if self.params.mmcif.hklout in (None, Auto, "auto"):
            if ("intensity.scale.value" in reflections) and (
                "intensity.scale.variance" in reflections
            ):
                filename = "scaled.cif"
                logger.info(
                    "Data appears to be scaled, setting mmcif.hklout = 'scaled.cif'"
                )
            else:
                filename = "integrated.cif"
                logger.info(
                    "Data appears to be unscaled, setting mmcif.hklout = 'integrated.cif'"
                )
        else:
            filename = self.params.mmcif.hklout

        # Add the block
        self._cif["dials"] = self.make_cif_block(experiments, reflections)

        loop_format_strings = {"_pdbx_diffrn_unmerged_refln": self._fmt}

        # Print to file
        if self.params.mmcif.compress and not filename.endswith(
            "." + self.params.mmcif.compress
        ):
            filename += "." + self.params.mmcif.compress
        if self.params.mmcif.compress == "gz":
            open_fn = gzip.open
        elif self.params.mmcif.compress == "bz2":
            open_fn = bz2.open
        elif self.params.mmcif.compress == "xz":
            open_fn = lzma.open
        else:
            open_fn = open
        with open_fn(filename, "wt") as fh:
            self._cif.show(out=fh, loop_format_strings=loop_format_strings)

        # Log
        logger.info("Wrote reflections to %s", filename)

    def make_cif_block(self, experiments, reflections):
        """Write the data to a cif block"""
        # Select reflections
        selection = reflections.get_flags(reflections.flags.integrated, all=True)
        reflections = reflections.select(selection)

        # Filter out bad variances and other issues, but don't filter on ice rings
        # or alter partialities.

        # Assumes you want to apply the lp and dqe corrections to sum and prf
        # Do we want to combine partials?
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
        cif_block["_audit.revision_id"] = 1
        cif_block["_audit.creation_method"] = dials_version
        cif_block["_audit.creation_date"] = datetime.date.today().isoformat()
        cif_block["_entry.id"] = "DIALS"
        # add software loop
        mmcif_software_header = (
            "_software.pdbx_ordinal",
            "_software.citation_id",
            "_software.name",  # as defined at [1]
            "_software.version",
            "_software.type",
            "_software.classification",
            "_software.description",
        )

        mmcif_citations_header = (
            "_citation.id",
            "_citation.journal_abbrev",
            "_citation.journal_volume",
            "_citation.journal_issue",
            "_citation.page_first",
            "_citation.page_last",
            "_citation.year",
            "_citation.title",
        )

        software_loop = iotbx.cif.model.loop(header=mmcif_software_header)
        citations_loop = iotbx.cif.model.loop(header=mmcif_citations_header)

        software_loop.add_row(
            (
                1,
                1,
                "DIALS",
                dials_version,
                "package",
                "data processing",
                "Data processing and integration within the DIALS software package",
            )
        )
        citations_loop.add_row(
            (
                1,
                "Acta Cryst. D",
                74,
                2,
                85,
                97,
                2018,
                "DIALS: implementation and evaluation of a new integration package",
            )
        )
        if "scale" in self.params.intensity:
            software_loop.add_row(
                (
                    2,
                    2,
                    "DIALS",
                    dials_version,
                    "program",
                    "data scaling",
                    "Data scaling and merging within the DIALS software package",
                )
            )
            citations_loop.add_row(
                (
                    2,
                    "Acta Cryst. D",
                    76,
                    4,
                    385,
                    399,
                    2020,
                    "Scaling diffraction data in the DIALS software package: algorithms and new approaches for multi-crystal scaling",
                )
            )
        cif_block.add_loop(software_loop)
        cif_block.add_loop(citations_loop)

        # Hard coding X-ray
        if self.params.mmcif.pdb_version == "v5_next":
            cif_block["_pdbx_diffrn_data_section.id"] = "dials"
            cif_block["_pdbx_diffrn_data_section.type_scattering"] = "x-ray"
            cif_block["_pdbx_diffrn_data_section.type_merged"] = "false"
            cif_block["_pdbx_diffrn_data_section.type_scaled"] = str(
                "scale" in self.params.intensity
            ).lower()

        # FIXME finish metadata addition - detector and source details needed
        # http://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Categories/index.html

        # Add source information;
        # _diffrn_source.pdbx_wavelength_list = (list of wavelengths)
        # _diffrn_source.source = (general class of source e.g. synchrotron)
        # _diffrn_source.type = (specific beamline or instrument e.g DIAMOND BEAMLINE I04)

        wls = []
        epochs = []
        for exp in experiments:
            wls.append(round(exp.beam.get_wavelength(), 5))
            epochs.append(exp.scan.get_epochs()[0])
        unique_wls = set(wls)
        cif_block["_exptl_crystal.id"] = 1  # links to crystal_id
        cif_block["_diffrn.id"] = 1  # links to diffrn_id
        cif_block["_diffrn.crystal_id"] = 1
        cif_block["_diffrn_source.diffrn_id"] = 1
        cif_block["_diffrn_source.pdbx_wavelength_list"] = ", ".join(
            str(w) for w in unique_wls
        )

        # Add detector information;
        # _diffrn_detector.detector  = (general class e.g. PIXEL, PLATE etc)
        # _diffrn_detector.pdbx_collection_date = (Date of collection yyyy-mm-dd)
        # _diffrn_detector.type = (full name of detector e.g. DECTRIS PILATUS3 2M)
        # One date is required, so if multiple just use the first date.
        min_epoch = min(epochs)
        date_str = time.strftime("%Y-%m-%d", time.gmtime(min_epoch))
        cif_block["_diffrn_detector.diffrn_id"] = 1
        cif_block["_diffrn_detector.pdbx_collection_date"] = date_str

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

        extra_items = {
            "scales": ("_pdbx_diffrn_unmerged_refln.scale_value", "%5.3f"),
            "intensity.scale.value": (
                "_pdbx_diffrn_unmerged_refln.intensity_meas",
                "%8.3f",
            ),
            "intensity.scale.sigma": (
                "_pdbx_diffrn_unmerged_refln.intensity_sigma",
                "%8.3f",
            ),
            "intensity.sum.value": (
                "_pdbx_diffrn_unmerged_refln.intensity_sum",
                "%8.3f",
            ),
            "intensity.sum.sigma": (
                "_pdbx_diffrn_unmerged_refln.intensity_sum_sigma",
                "%8.3f",
            ),
            "intensity.prf.value": (
                "_pdbx_diffrn_unmerged_refln.intensity_prf",
                "%8.3f",
            ),
            "intensity.prf.sigma": (
                "_pdbx_diffrn_unmerged_refln.intensity_prf_sigma",
                "%8.3f",
            ),
            "angle": ("_pdbx_diffrn_unmerged_refln.scan_angle_reflection", "%7.4f"),
            "partiality": ("_pdbx_diffrn_unmerged_refln.partiality", "%7.4f"),
        }

        variables_present = []
        if "scale" in self.params.intensity:
            reflections["scales"] = 1.0 / reflections["inverse_scale_factor"]
            reflections["intensity.scale.sigma"] = flex.sqrt(
                reflections["intensity.scale.variance"]
            )
            variables_present.extend(
                ["scales", "intensity.scale.value", "intensity.scale.sigma"]
            )
        if "sum" in self.params.intensity:
            reflections["intensity.sum.sigma"] = flex.sqrt(
                reflections["intensity.sum.variance"]
            )
            variables_present.extend(["intensity.sum.value", "intensity.sum.sigma"])
        if "profile" in self.params.intensity:
            reflections["intensity.prf.sigma"] = flex.sqrt(
                reflections["intensity.prf.variance"]
            )
            variables_present.extend(["intensity.prf.value", "intensity.prf.sigma"])

        # Should always exist
        reflections["angle"] = reflections["xyzcal.mm"].parts()[2] * RAD2DEG
        variables_present.extend(["angle"])

        if "partiality" in reflections:
            variables_present.extend(["partiality"])

        for name in variables_present:
            if name in reflections:
                header += (extra_items[name][0],)
                self._fmt += " " + extra_items[name][1]

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
                assert_is_not_unique_set_under_symmetry=False,
            )
            merged_block = iotbx.cif.model.block()
            merged_block["_reflns.pdbx_ordinal"] = 1
            merged_block["_reflns.pdbx_diffrn_id"] = 1
            merged_block["_reflns.entry_id"] = "DIALS"
            merged_data = result.as_cif_block()
            merged_block.update(merged_data)
            cif_block.update(merged_block)

        # Write the crystal information
        # if v5, that's all so return
        if self.params.mmcif.pdb_version == "v5":
            return cif_block
        # continue if v5_next
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

        expid_to_scan_id = {exp.identifier: i + 1 for i, exp in enumerate(experiments)}

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

        _, _, _, _, z0, z1 = reflections["bbox"].parts()
        h, k, l = [
            hkl.iround() for hkl in reflections["miller_index"].as_vec3_double().parts()
        ]
        # make scan id consistent with header as defined above
        scan_id = flex.int(reflections.size(), 0)
        for id_ in reflections.experiment_identifiers().keys():
            expid = reflections.experiment_identifiers()[id_]
            sel = reflections["id"] == id_
            scan_id.set_selected(sel, expid_to_scan_id[expid])

        loop_values = [
            flex.size_t_range(1, len(reflections) + 1),
            scan_id,
            z0,
            z1,
            h,
            k,
            l,
        ] + [reflections[name] for name in variables_present]
        cif_loop = iotbx.cif.model.loop(data=dict(zip(header, loop_values)))
        cif_block.add_loop(cif_loop)

        return cif_block

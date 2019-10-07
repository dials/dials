# coding: utf-8

from __future__ import absolute_import, division, print_function

import logging
import time
from collections import OrderedDict, Counter

import dials.util.ext
from dials.array_family import flex
from dials.util.version import dials_version
from dials.util.filter_reflections import filter_reflection_table
from dials.util.batch_handling import (
    calculate_batch_offsets,
    assign_batches_to_reflections,
    get_image_ranges,
)
from dials.util.multi_dataset_handling import (
    assign_unique_identifiers,
    parse_multiple_datasets,
)
from iotbx import mtz
from scitbx import matrix
from libtbx import env, Auto

try:
    from math import isclose
except ImportError:
    # Python 3 backport
    def isclose(a, b, rel_tol=1e-09, abs_tol=0.0):
        return abs(a - b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)


logger = logging.getLogger(__name__)


class MTZWriterBase(object):
    """Helper for adding metadata, crystals and datasets to an mtz file object."""

    def __init__(self, space_group, unit_cell=None):
        """If a unit cell is provided, will be used as default unless specified
        for each crystal."""
        mtz_file = mtz.object()
        mtz_file.set_title("From %s" % env.dispatcher_name)
        date_str = time.strftime("%Y-%m-%d at %H:%M:%S %Z")
        if time.strftime("%Z") != "GMT":
            date_str += time.strftime("  (%Y-%m-%d at %H:%M:%S %Z)", time.gmtime())
        mtz_file.add_history("From %s, run on %s" % (dials_version(), date_str))
        mtz_file.set_space_group_info(space_group.info())
        self.mtz_file = mtz_file
        if unit_cell:
            self.unit_cell = unit_cell
        self.current_crystal = None
        self.current_dataset = None
        self.n_crystals = 0
        self.n_datasets = 0

    def add_crystal(self, crystal_name=None, project_name=None, unit_cell=None):
        """Add a crystal to the mtz file object."""
        if not unit_cell:
            if not self.unit_cell:
                raise ValueError("Unit cell must be provided.")
            else:
                unit_cell = self.unit_cell
        if not crystal_name:
            crystal_name = "crystal_%s" % str(self.n_crystals + 1)
        if not project_name:
            project_name = "DIALS"
        self.current_crystal = self.mtz_file.add_crystal(
            crystal_name, project_name, unit_cell.parameters()
        )
        self.n_crystals += 1

    def add_empty_dataset(self, wavelength, name=None):
        """Add an empty dataset object to the mtz file."""
        if not name:
            name = "FROMDIALS"
        self.current_dataset = self.current_crystal.add_dataset(name, wavelength)
        self.n_datasets += 1


class MergedMTZWriter(MTZWriterBase):
    """Mtz writer for merged data."""

    def add_dataset(
        self,
        merged_array,
        anom_array=None,
        amplitudes=None,
        anom_amplitudes=None,
        suffix=None,
    ):
        """Add merged data to the most recent dataset.

        Args:
            merged_array: A merged miller array of IMEAN intensities
            wavelength: The wavelength of the dataset
            anom_array (Optional): An anomalous merged miller array
            amplitudes (Optional): A merged miller array of amplitudes
            anom_amplitudes (Optional): An anomalous merged array of amplitudes
            suffix (Optional[str]): Column name suffix to use for this dataset.
        """
        if not suffix:
            suffix = ""
        self.current_dataset.add_miller_array(merged_array, "IMEAN" + suffix)
        if anom_array:
            self.current_dataset.add_miller_array(anom_array, "I" + suffix)
            self.current_dataset.add_miller_array(
                anom_array.multiplicities(),
                column_root_label="N" + suffix,
                column_types="I",
            )
        if amplitudes:
            self.current_dataset.add_miller_array(amplitudes, "F" + suffix)
        if anom_amplitudes:
            self.current_dataset.add_miller_array(anom_amplitudes, "F" + suffix)


class MADMergedMTZWriter(MergedMTZWriter):
    """Mtz writer for multi-wavelength merged data."""

    def add_dataset(
        self,
        merged_array,
        anom_array=None,
        amplitudes=None,
        anom_amplitudes=None,
        suffix=None,
    ):
        if not suffix:
            suffix = "_WAVE%s" % str(self.n_datasets)
        super(MADMergedMTZWriter, self).add_dataset(
            merged_array, anom_array, amplitudes, anom_amplitudes, suffix
        )


class UnmergedMTZWriter(MTZWriterBase):
    def add_batch_list(
        self,
        image_range,
        experiment,
        wavelength,
        dataset_id,
        batch_offset,
        force_static_model,
    ):
        """Add batch metadata to the mtz file."""

        # Recalculate useful numbers and references here
        n_batches = image_range[1] - image_range[0] + 1
        phi_start = flex.float(n_batches, 0)
        phi_range = flex.float(n_batches, 0)
        umat_array = flex.float(flex.grid(n_batches, 9))
        cell_array = flex.float(flex.grid(n_batches, 6))

        U = matrix.sqr(experiment.crystal.get_U())
        if experiment.goniometer is not None:
            F = matrix.sqr(experiment.goniometer.get_fixed_rotation())
        else:
            F = matrix.sqr((1, 0, 0, 0, 1, 0, 0, 0, 1))

        i0 = image_range[0]
        for i in range(n_batches):
            if experiment.scan:
                phi_start[i], phi_range[i] = experiment.scan.get_image_oscillation(
                    i + i0
                )

            # unit cell (this is fine) and the what-was-refined-flags hardcoded
            # take time-varying parameters from the *end of the frame* unlikely to
            # be much different at the end - however only exist if scan-varying
            # refinement was used
            if not force_static_model and experiment.crystal.num_scan_points > 0:
                # Get the index of the image in the sequence e.g. first => 0, second => 1
                image_index = i + i0 - experiment.scan.get_image_range()[0]
                _unit_cell = experiment.crystal.get_unit_cell_at_scan_point(image_index)
                _U = matrix.sqr(experiment.crystal.get_U_at_scan_point(image_index))
            else:
                _unit_cell = experiment.crystal.get_unit_cell()
                _U = U

            # apply the fixed rotation to this to unify matrix definitions - F * U
            # was what was used in the actual prediction: U appears to be stored
            # as the transpose?! At least is for Mosflm...
            #
            # FIXME Do we need to apply the setting rotation here somehow? i.e. we have
            # the U.B. matrix assuming that the axis is equal to S * axis_datum but
            # here we are just giving the effective axis so at scan angle 0 this will
            # not be correct... FIXME 2 not even sure we can express the stack of
            # matrices S * R * F * U * B in MTZ format?... see [=A=] below
            _U = matrix.sqr(dials.util.ext.dials_u_to_mosflm(F * _U, _unit_cell))

            # FIXME need to get what was refined and what was constrained from the
            # crystal model - see https://github.com/dials/dials/issues/355
            _unit_cell_params = _unit_cell.parameters()
            for j in range(6):
                cell_array[i, j] = _unit_cell_params[j]
            _U_t_elements = _U.transpose().elems
            for j in range(9):
                umat_array[i, j] = _U_t_elements[j]

        # We ignore panels beyond the first one, at the moment
        panel = experiment.detector[0]
        panel_size = panel.get_image_size()
        panel_distance = panel.get_directed_distance()

        if experiment.goniometer:
            axis = flex.float(experiment.goniometer.get_rotation_axis())
        else:
            axis = flex.float((0.0, 0.0, 0.0))

        # FIXME hard-coded assumption on idealized beam vector below... this may be
        # broken when we come to process data from a non-imgCIF frame
        s0n = flex.float(matrix.col(experiment.beam.get_s0()).normalize().elems)

        # get the mosaic spread though today it may not actually be set - should
        # this be in the BATCH headers?
        try:
            mosaic = experiment.crystal.get_mosaicity()
        except AttributeError:
            mosaic = 0.0

        # Jump into C++ to do the rest of the work
        dials.util.ext.add_dials_batches(
            self.mtz_file,
            dataset_id,
            image_range,
            batch_offset,
            wavelength,
            mosaic,
            phi_start,
            phi_range,
            cell_array,
            umat_array,
            panel_size,
            panel_distance,
            axis,
            s0n,
        )

    def write_columns(self, integrated_data):
        """Write the column definitions AND data to the current dataset."""

        # now create the actual data structures - first keep a track of the columns

        # H K L M/ISYM BATCH I SIGI IPR SIGIPR FRACTIONCALC XDET YDET ROT WIDTH
        # LP MPART FLAG BGPKRATIOS

        # gather the required information for the reflection file

        nref = len(integrated_data["miller_index"])
        assert nref
        xdet, ydet, _ = [
            flex.double(x) for x in integrated_data["xyzobs.px.value"].parts()
        ]

        # now add column information...

        # FIXME add DIALS_FLAG which can include e.g. was partial etc.

        type_table = {
            "H": "H",
            "K": "H",
            "L": "H",
            "I": "J",
            "SIGI": "Q",
            "IPR": "J",
            "SIGIPR": "Q",
            "BG": "R",
            "SIGBG": "R",
            "XDET": "R",
            "YDET": "R",
            "BATCH": "B",
            "BGPKRATIOS": "R",
            "WIDTH": "R",
            "MPART": "I",
            "M_ISYM": "Y",
            "FLAG": "I",
            "LP": "R",
            "FRACTIONCALC": "R",
            "ROT": "R",
            "QE": "R",
        }

        # derive index columns from original indices with
        #
        # from m.replace_original_index_miller_indices
        #
        # so all that is needed now is to make space for the reflections - fill with
        # zeros...

        self.mtz_file.adjust_column_array_sizes(nref)
        self.mtz_file.set_n_reflections(nref)
        dataset = self.current_dataset

        # assign H, K, L, M_ISYM space
        for column in "H", "K", "L", "M_ISYM":
            dataset.add_column(column, type_table[column]).set_values(
                flex.double(nref, 0.0).as_float()
            )

        self.mtz_file.replace_original_index_miller_indices(
            integrated_data["miller_index_rebase"]
        )

        dataset.add_column("BATCH", type_table["BATCH"]).set_values(
            integrated_data["batch"].as_double().as_float()
        )

        # if intensity values used in scaling exist, then just export these as I, SIGI
        if "intensity.scale.value" in integrated_data:
            I_scaling = integrated_data["intensity.scale.value"]
            V_scaling = integrated_data["intensity.scale.variance"]
            # Trap negative variances
            assert V_scaling.all_gt(0)
            dataset.add_column("I", type_table["I"]).set_values(I_scaling.as_float())
            dataset.add_column("SIGI", type_table["SIGI"]).set_values(
                flex.sqrt(V_scaling).as_float()
            )
            dataset.add_column("SCALEUSED", "R").set_values(
                integrated_data["inverse_scale_factor"].as_float()
            )
            dataset.add_column("SIGSCALEUSED", "R").set_values(
                flex.sqrt(integrated_data["inverse_scale_factor_variance"]).as_float()
            )
        else:
            if "intensity.prf.value" in integrated_data:
                if "intensity.sum.value" in integrated_data:
                    col_names = ("IPR", "SIGIPR")
                else:
                    col_names = ("I", "SIGI")
                I_profile = integrated_data["intensity.prf.value"]
                V_profile = integrated_data["intensity.prf.variance"]
                # Trap negative variances
                assert V_profile.all_gt(0)
                dataset.add_column(col_names[0], type_table["I"]).set_values(
                    I_profile.as_float()
                )
                dataset.add_column(col_names[1], type_table["SIGI"]).set_values(
                    flex.sqrt(V_profile).as_float()
                )
            if "intensity.sum.value" in integrated_data:
                I_sum = integrated_data["intensity.sum.value"]
                V_sum = integrated_data["intensity.sum.variance"]
                # Trap negative variances
                assert V_sum.all_gt(0)
                dataset.add_column("I", type_table["I"]).set_values(I_sum.as_float())
                dataset.add_column("SIGI", type_table["SIGI"]).set_values(
                    flex.sqrt(V_sum).as_float()
                )
        if (
            "background.sum.value" in integrated_data
            and "background.sum.variance" in integrated_data
        ):
            bg = integrated_data["background.sum.value"]
            varbg = integrated_data["background.sum.variance"]
            assert (varbg >= 0).count(False) == 0
            sigbg = flex.sqrt(varbg)
            dataset.add_column("BG", type_table["BG"]).set_values(bg.as_float())
            dataset.add_column("SIGBG", type_table["SIGBG"]).set_values(
                sigbg.as_float()
            )

        dataset.add_column("FRACTIONCALC", type_table["FRACTIONCALC"]).set_values(
            integrated_data["fractioncalc"].as_float()
        )

        dataset.add_column("XDET", type_table["XDET"]).set_values(xdet.as_float())
        dataset.add_column("YDET", type_table["YDET"]).set_values(ydet.as_float())
        dataset.add_column("ROT", type_table["ROT"]).set_values(
            integrated_data["ROT"].as_float()
        )
        if "lp" in integrated_data:
            dataset.add_column("LP", type_table["LP"]).set_values(
                integrated_data["lp"].as_float()
            )
        if "qe" in integrated_data:
            dataset.add_column("QE", type_table["QE"]).set_values(
                integrated_data["qe"].as_float()
            )
        elif "dqe" in integrated_data:
            dataset.add_column("QE", type_table["QE"]).set_values(
                integrated_data["dqe"].as_float()
            )
        else:
            dataset.add_column("QE", type_table["QE"]).set_values(
                flex.double(nref, 1.0).as_float()
            )


def export_mtz(integrated_data, experiment_list, params):
    """Export data from integrated_data corresponding to experiment_list to an
    MTZ file hklout."""

    # if mtz filename is auto, then choose scaled.mtz or integrated.mtz
    if params.mtz.hklout in (None, Auto, "auto"):
        if ("intensity.scale.value" in integrated_data) and (
            "intensity.scale.variance" in integrated_data
        ):
            params.mtz.hklout = "scaled.mtz"
            logger.info(
                "Data appears to be scaled, setting mtz.hklout = 'scaled_unmerged.mtz'"
            )
        else:
            params.mtz.hklout = "integrated.mtz"
            logger.info(
                "Data appears to be unscaled, setting mtz.hklout = 'integrated.mtz'"
            )

    # First get the experiment identifier information out of the data
    expids_in_table = integrated_data.experiment_identifiers()
    if not list(expids_in_table.keys()):
        reflection_tables = parse_multiple_datasets([integrated_data])
        experiment_list, refl_list = assign_unique_identifiers(
            experiment_list, reflection_tables
        )
        integrated_data = flex.reflection_table()
        for reflections in refl_list:
            integrated_data.extend(reflections)
        expids_in_table = integrated_data.experiment_identifiers()
    integrated_data.assert_experiment_identifiers_are_consistent(experiment_list)
    expids_in_list = list(experiment_list.identifiers())

    # Convert experiment_list to a real python list or else identity assumptions
    # fail like:
    #   assert experiment_list[0] is experiment_list[0]
    # And assumptions about added attributes break
    experiment_list = list(experiment_list)

    # Validate multi-experiment assumptions
    if len(experiment_list) > 1:
        # All experiments should match crystals, or else we need multiple crystals/datasets
        if not all(
            x.crystal == experiment_list[0].crystal for x in experiment_list[1:]
        ):
            logger.warning(
                "Warning: Experiment crystals differ. Using first experiment crystal for file-level data."
            )

        wavelengths = match_wavelengths(experiment_list)
        if len(wavelengths.keys()) > 1:
            logger.info(
                "Multiple wavelengths found: \n%s",
                "\n".join(
                    "  Wavlength: %.5f, experiment numbers: %s "
                    % (k, ",".join(map(str, v)))
                    for k, v in wavelengths.items()
                ),
            )
    else:
        wavelengths = OrderedDict({experiment_list[0].beam.get_wavelength(): [0]})

    # also only work correctly with one panel (for the moment)
    if any(len(experiment.detector) != 1 for experiment in experiment_list):
        logger.warning("Warning: Ignoring multiple panels in output MTZ")

    # Clean up the data with the passed in options
    integrated_data = filter_reflection_table(
        integrated_data,
        intensity_choice=params.intensity,
        partiality_threshold=params.mtz.partiality_threshold,
        combine_partials=params.mtz.combine_partials,
        min_isigi=params.mtz.min_isigi,
        filter_ice_rings=params.mtz.filter_ice_rings,
        d_min=params.mtz.d_min,
    )

    # get batch offsets and image ranges - even for scanless experiments
    batch_offsets = [
        expt.scan.get_batch_offset()
        for expt in experiment_list
        if expt.scan is not None
    ]
    unique_offsets = set(batch_offsets)
    if len(set(unique_offsets)) <= 1:
        logger.debug("Calculating new batches")
        batch_offsets = calculate_batch_offsets(experiment_list)
        batch_starts = [
            e.scan.get_image_range()[0] if e.scan else 0 for e in experiment_list
        ]
        effective_offsets = [o + s for o, s in zip(batch_offsets, batch_starts)]
        unique_offsets = set(effective_offsets)
    else:
        logger.debug("Keeping existing batches")
    image_ranges = get_image_ranges(experiment_list)
    if len(unique_offsets) != len(batch_offsets):

        raise ValueError(
            "Duplicate batch offsets detected: %s"
            % ", ".join(
                str(item) for item, count in Counter(batch_offsets).items() if count > 1
            )
        )

    # Create the mtz file
    mtz_writer = UnmergedMTZWriter(experiment_list[0].crystal.get_space_group())

    # FIXME TODO for more than one experiment into an MTZ file:
    #
    # - add an epoch (or recover an epoch) from the scan and add this as an extra
    #   column to the MTZ file for scaling, so we know that the two lattices were
    #   integrated at the same time
    # ✓ decide a sensible BATCH increment to apply to the BATCH value between
    #   experiments and add this

    for id_ in expids_in_table.keys():
        # Grab our subset of the data
        loc = expids_in_list.index(
            expids_in_table[id_]
        )  # get strid and use to find loc in list
        experiment = experiment_list[loc]
        if len(list(wavelengths.keys())) > 1:
            for i, (wl, exps) in enumerate(wavelengths.items()):
                if loc in exps:
                    wavelength = wl
                    dataset_id = i + 1
                    break
        else:
            wavelength = list(wavelengths.keys())[0]
            dataset_id = 1
        reflections = integrated_data.select(integrated_data["id"] == id_)
        batch_offset = batch_offsets[loc]
        image_range = image_ranges[loc]
        reflections = assign_batches_to_reflections([reflections], [batch_offset])[0]
        experiment.data = dict(reflections)

        # Do any crystal transformations for the experiment
        cb_op_to_ref = (
            experiment.crystal.get_space_group()
            .info()
            .change_of_basis_op_to_reference_setting()
        )
        experiment.crystal = experiment.crystal.change_basis(cb_op_to_ref)
        experiment.data["miller_index_rebase"] = cb_op_to_ref.apply(
            experiment.data["miller_index"]
        )

        s0n = matrix.col(experiment.beam.get_s0()).normalize().elems
        logger.debug("Beam vector: %.4f %.4f %.4f" % s0n)

        mtz_writer.add_batch_list(
            image_range,
            experiment,
            wavelength,
            dataset_id,
            batch_offset=batch_offset,
            force_static_model=params.mtz.force_static_model,
        )

        # Create the batch offset array. This gives us an experiment (id)-dependent
        # batch offset to calculate the correct batch from image number.
        experiment.data["batch_offset"] = flex.int(
            len(experiment.data["id"]), batch_offset
        )

        # Calculate whether we have a ROT value for this experiment, and set the column
        _, _, z = experiment.data["xyzcal.px"].parts()
        if experiment.scan:
            experiment.data["ROT"] = experiment.scan.get_angle_from_array_index(z)
        else:
            experiment.data["ROT"] = z

    mtz_writer.add_crystal(
        params.mtz.crystal_name, unit_cell=experiment_list[0].crystal.get_unit_cell()
    )  # Note: add unit cell here as may have changed basis since creating mtz.
    # For multi-wave unmerged mtz, we add an empty dataset for each wavelength,
    # but only write the data into the final dataset (for unmerged the batches
    # link the unmerged data to the individual wavelengths).
    for wavelength in wavelengths:
        mtz_writer.add_empty_dataset(wavelength)

    # Combine all of the experiment data columns before writing
    combined_data = {k: v.deep_copy() for k, v in experiment_list[0].data.items()}
    for experiment in experiment_list[1:]:
        for k, v in experiment.data.items():
            combined_data[k].extend(v)
    # ALL columns must be the same length
    assert len({len(v) for v in combined_data.values()}) == 1, "Column length mismatch"
    assert len(combined_data["id"]) == len(
        integrated_data["id"]
    ), "Lost rows in split/combine"

    # Write all the data and columns to the mtz file
    mtz_writer.write_columns(combined_data)

    logger.info(
        "Saving {} integrated reflections to {}".format(
            len(combined_data["id"]), params.mtz.hklout
        )
    )
    mtz_file = mtz_writer.mtz_file
    mtz_file.write(params.mtz.hklout)

    return mtz_file


def match_wavelengths(experiments):
    """Create a dictionary matching wavelength to experiments (index in list)"""
    wavelengths = OrderedDict()
    for i, x in enumerate(experiments):
        w = x.beam.get_wavelength()
        matches = [isclose(w, k, rel_tol=1e-4) for k in wavelengths]
        if not any(matches):
            wavelengths[w] = [i]
        else:
            match_w = list(wavelengths.keys())[matches.index(True)]
            wavelengths[match_w].append(i)
    return wavelengths

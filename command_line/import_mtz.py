"""
dials.import_mtz, a tool for creating DIALS working data formats (reflection
tables and experiments) from an MTZ file.
"""

from __future__ import annotations

import logging
import math
import sys
from collections import OrderedDict
from dataclasses import dataclass
from typing import List, Tuple

import numpy as np
from scipy.spatial.transform import Rotation
from scipy.stats import linregress

import cctbx
import iotbx.mtz
import libtbx.phil
import scitbx.array_family
from dxtbx import flumpy
from dxtbx.model import Beam, Crystal, Detector, Experiment, ExperimentList, Scan
from dxtbx.model.goniometer import GoniometerFactory
from rstbx.cftbx.coordinate_frame_helpers import align_reference_frame
from scitbx import matrix

from dials.algorithms.spot_prediction import ray_intersection
from dials.array_family import flex
from dials.util import log, show_mail_handle_errors, tabulate
from dials.util.exclude_images import set_initial_valid_image_ranges
from dials.util.multi_dataset_handling import generate_experiment_identifiers
from dials.util.options import OptionParser

logger = logging.getLogger("dials.command_line.import_mtz")

help_message = """

This program makes DIALS datafiles from an MTZ input.

More specifically, it is intended to allow importing of integrated MTZ files
for data reduction in DIALS or dependent software, e.g. xia2.multiplex.
As not all information can be contained within mtz files, there may be some
less common program options that are not possible with these files imported
from mtz.

Ideally, the image template can be provided with image.template=, which allows
reading of the (unrefined) detector models, goniometer, beam and scan, and the
crystal model can be determined from the MTZ data.
If the template is not provided, then the models are constructed from the header
data. It may be necessary to provide the scan oscillation step.


"""

phil_scope = libtbx.phil.parse(
    """

  input {

    use_datasets = None
      .type = ints
      .help = "Select a subset of datasets from the input MTZ file."

    oscillation_step = None
      .type = float
      .help = "The rotational width of an image, used to convert rotation"
              "angle to image number."

    labels {
      intensity = "I,SIGI"
        .type = str
        .help = "column labels (Iname,SIGIname) for the intensity array"
      xpos = "XDET"
        .type = str
        .help = "column label for spot x-coordinate on the detector"
      ypos = "YDET"
        .type = str
        .help = "column label for spot y-coordinate on the detector"
      phi = "ROT"
        .type = str
        .help = "column label for rotation angle."
      partiality = "FRACTIONCALC"
        .type = str
        .help = "column label for partiality."

    }

    images {
        template = None
            .type = str
            .help = "The image sequence template"
        image_range = None
            .type = ints(size=2)
            .help = "First and last image in the template"
    }

    detector {
        pixel_size = None
            .type = float(value_min=0.0)
            .help = "Pixel size in mm"
    }

    include scope dxtbx.model.goniometer.goniometer_phil_scope

  }

  output {

    reflections = imported_mtz.refl
      .type = str
      .help = "Name of the output reflections file"

    experiments = imported_mtz.expt
      .type = str
      .help = "Name of the output experiments file"

  }
""",
    process_includes=True,
)


def detector_info_from_batch_headers(
    unmerged_mtz: iotbx.mtz.object,
) -> Tuple[Tuple[int, int], float]:
    b0 = unmerged_mtz.batches()[0]
    panel_size = int(b0.detlm()[1]), int(b0.detlm()[3])
    panel_distance = b0.dx()[0]
    return panel_size, panel_distance


def mosflm_B_matrix(unit_cell: cctbx.uctbx.unit_cell) -> scitbx.matrix:
    rp = unit_cell.reciprocal_parameters()
    p = unit_cell.parameters()
    d2r = math.pi / 180.0
    B = matrix.sqr(
        (
            rp[0],
            rp[1] * math.cos(d2r * rp[5]),
            rp[2] * math.cos(d2r * rp[4]),
            0,
            rp[1] * math.sin(d2r * rp[5]),
            -rp[2] * math.sin(d2r * rp[4]) * math.cos(d2r * p[3]),
            0,
            0,
            1.0 / p[2],
        )
    )
    return B


def scan_info_from_batch_headers(
    unmerged_mtz: iotbx.mtz.object, verbose: bool = True
) -> OrderedDict:
    batches = unmerged_mtz.batches()

    # Determine rotation matrix to convert to the DIALS frame
    scanax = matrix.col(batches[0].scanax())
    if scanax.length() == 0.0:
        # Default for MOSFLM, which does not set scanax
        scanax = matrix.col((0, 0, 1))
    M = align_reference_frame(
        matrix.col(batches[0].source()),
        matrix.col((0, 0, 1)),
        matrix.col(scanax),
        matrix.col((1, 0, 0)),
    )

    scans = OrderedDict(
        {
            1: {
                "start_image": 1,
                "end_image": None,
                "batch_begin": batches[0].num(),
                "batch_end": None,
                "angle_begin": batches[0].phistt(),
                "angle_end": None,
                "umat_start": batches[0].umat(),
                "s0n": batches[0].so(),
                "transform_to_dials": M,
            }
        }
    )

    scan_no = 1
    phi_end = batches[0].phiend()
    last_batch = batches[0].num()

    for b in batches[1:]:
        if abs(b.phistt() - phi_end) > 0.0001:

            # Determine rotation matrix to convert to the DIALS frame
            M = align_reference_frame(
                matrix.col(b[0].source()),
                matrix.col((0, 0, -1)),
                matrix.col(b[0].scanax()),
                matrix.col((1, 0, 0)),
            )

            scans[scan_no]["angle_end"] = phi_end
            scans[scan_no]["batch_end"] = last_batch
            scans[scan_no]["end_image"] = last_batch - scans[scan_no]["batch_begin"] + 1

            scan_no += 1
            scans[scan_no] = {
                "start_image": 1,
                "end_image": None,
                "batch_begin": b.num(),
                "batch_end": None,
                "angle_begin": b.phistt(),
                "angle_end": None,
                "umat_start": batches[0].umat(),
                "s0n": batches[0].so(),
                "transform_to_dials": M,
            }

        phi_end = b.phiend()
        last_batch = b.num()

    scans[scan_no]["angle_end"] = phi_end
    scans[scan_no]["batch_end"] = last_batch
    scans[scan_no]["end_image"] = last_batch - scans[scan_no]["batch_begin"] + 1

    if verbose:
        header = ["Scan number", "image range", "batch range", "scan range"]
        rows = []
        for n, s in scans.items():
            ifirst, ilast = s["start_image"], s["end_image"]
            bfirst, blast = s["batch_begin"], s["batch_end"]
            afirst, alast = s["angle_begin"], s["angle_end"]
            rows.append(
                [
                    f"{n}",
                    f"{ifirst} : {ilast}",
                    f"{bfirst} : {blast}",
                    f"{afirst} : {alast}",
                ]
            )
        logger.info("Scan information determined from batch headers:")
        logger.info(tabulate(rows, header))

    return scans


@dataclass
class MTZColumnData:
    intensities: cctbx.miller.array = None
    profile_intensities: cctbx.miller.array = None
    batches: scitbx.array_family.flex.int = None
    indices: cctbx.array_family.flex.miller_index = None
    xpos: scitbx.array_family.flex.double = None
    ypos: scitbx.array_family.flex.double = None
    phi: scitbx.array_family.flex.double = None
    scale_used: scitbx.array_family.flex.double = None
    partiality: scitbx.array_family.flex.double = None
    bg: scitbx.array_family.flex.double = None
    bgsig: scitbx.array_family.flex.double = None
    lp: scitbx.array_family.flex.double = None
    qe: scitbx.array_family.flex.double = None


def extract_data_from_mtz(
    unmerged_mtz: iotbx.mtz.object, params: libtbx.phil.scope_extract
) -> MTZColumnData:

    miller_arrays = unmerged_mtz.as_miller_arrays(
        merge_equivalents=False, anomalous=True
    )
    labels = [label for array in miller_arrays for label in array.info().labels]
    logger.info(f"\nData arrays found in mtz file: {', '.join(labels)}")

    data = MTZColumnData()
    data.indices = unmerged_mtz.extract_original_index_miller_indices()
    for array in miller_arrays:
        labels = array.info().labels
        if labels == [params.input.labels.xpos]:
            data.xpos = array.data()
        elif labels == [params.input.labels.ypos]:
            data.ypos = array.data()
        elif labels == [params.input.labels.phi]:
            data.phi = array.data()  # this is xyzcal.px
        elif labels == [params.input.labels.partiality]:
            data.partiality = array.data()
            # want to transform from rotation angle to image number
        elif labels == params.input.labels.intensity.split(","):
            data.intensities = array
        elif labels == ["IPR", "SIGIPR"]:
            data.profile_intensities = array
        elif labels == ["BATCH"]:
            data.batches = array.data()
        elif labels == ["BG"]:
            data.bg = array.data()
        elif labels == ["SIGBG"]:
            data.bgsig = array.data()
        elif labels == ["LP"]:
            data.lp = array.data()
        elif labels == ["QE"]:
            data.qe = array.data()
        elif labels == ["SCALEUSED"]:
            data.scale_used = array.data()

    if not data.intensities:
        raise ValueError(
            f"Intensities not found using labels {params.input.labels.intensity}"
            + "\nIntensity column labels can be specified with labels.intensity='Iname,SIGIname'"
        )
    if not data.batches:
        raise ValueError("Batch values not found")
    if not data.phi:
        raise ValueError(
            f"""Phi/Rot values not found using labels {params.input.labels.phi}
Phi/Rot column label can be specified with labels.phi='PHIname'"""
        )
    if data.batches.size() != data.intensities.data().size():
        raise ValueError("Batch and intensity array sizes do not match")
    return data


def create_reflection_table_from_mtz(
    unmerged_mtz: iotbx.mtz.object,
    params: libtbx.phil.scope_extract,
    experiments: ExperimentList,
    batch_offset: int = None,
) -> flex.reflection_table:
    """Extract reflection table data"""

    # xds integrate.hkl file: IOBS,SIGMA are profile fitted values. Intensity
    # is LP corrected assuming unpolarised beam. Polarisation correction carried
    # out in correct step.
    # xds integrate.hkl after converting with pointless: I,SIGI are profile fitted
    # presume LP correction not applied.

    # Dials integrated.mtz - default I/SIGI = summation IPR,SIGIPR - profile fitted
    # this is what aimless expects in an input mtz. but what about a dials integrated
    # dataset with only summation intensities - looks same as xds profile only?

    # Dials scaled.mtz - just I/SIGI are the scaled intensities i.e. ready for merging.
    # uniqueness is presence of SCALEUSED column.

    data = extract_data_from_mtz(unmerged_mtz, params)
    scans = scan_info_from_batch_headers(unmerged_mtz, verbose=False)

    if len(scans) > 1:
        raise ValueError(
            "Only importing MTZ files containing one sweep is currently supported"
        )

    # The reflection data
    table = flex.reflection_table()
    table["miller_index"] = data.indices
    table["d"] = data.intensities.d_spacings().data()
    table["partiality"] = (
        data.partiality if data.partiality else flex.double(table.size(), 1.0)
    )

    # extract the intensities. In dials, we work on raw intensities.
    # In mtz files, the intensities already have corrections applied. i.e
    # if we have IPR and LP, then intensity.prf.value = IPR/LP

    # First, if the SCALEUSED column is there, indicates scaled data
    # In that case, I,SIGI are scaled intensities
    if data.scale_used:
        logger.info(
            """\nThe intensity data with labels I,SIGI are interpeted as being
scaled intensity data. Any further scaling correction will be made in
addition to the existing correction (SCALEUSED)."""
        )
        # Assertion is that we want to use the scaled intensities.
        # Set the scaled intensities as 'intensity.prf.value'
        table["intensity.prf.value"] = data.intensities.data()
        table["intensity.prf.variance"] = flex.pow2(data.intensities.sigmas())
        table["scaleused"] = data.scale_used
        # now these have had corrections like lp applied already. So won't
        # be wanting to apply these again, therefore don't set these columns.
        table.set_flags(flex.bool(table.size(), True), table.flags.integrated_prf)

    # FIXME check these assumptions for aimless, xds, etc data
    # now see whether two or one sets of intensities were found
    elif data.profile_intensities and data.intensities:
        logger.info(
            """\nThe intensity data with labels I,SIGI are interpreted as being
summation integrated intensities.
The data with labels IPR,SIGIPR are interpreted as being
profile-fitted integrated intensities."""
        )
        table["intensity.sum.value"] = data.intensities.data()
        table["intensity.sum.variance"] = flex.pow2(data.intensities.sigmas())
        table["intensity.prf.value"] = data.profile_intensities.data()
        table["intensity.prf.variance"] = flex.pow2(data.profile_intensities.sigmas())
        table.set_flags(flex.bool(table.size(), True), table.flags.integrated_sum)
        table.set_flags(flex.bool(table.size(), True), table.flags.integrated_prf)
        ## In DIALS, we work on raw values, so need to reverse the lp, qe, partiality corrections
        correction = flex.double(table.size(), 1.0)
        if data.lp:
            table["lp"] = data.lp
            assert data.lp.all_gt(0.0)
            correction /= data.lp
        if data.qe:
            table["qe"] = data.qe
            correction *= data.qe
        if data.qe or data.lp:
            table["intensity.sum.value"] *= correction
            table["intensity.sum.variance"] *= correction**2
            table["intensity.prf.value"] *= correction
            table["intensity.prf.variance"] *= correction**2
        if data.partiality:
            table["intensity.sum.value"] *= data.partiality
            table["intensity.sum.variance"] *= data.partiality**2

    else:  # only one intensity set, assume profile intensities.
        logger.info(
            """The data with labels I,SIGI are interpreted as being
profile-fitted integrated intensities."""
        )
        table["intensity.prf.value"] = data.intensities.data()
        table["intensity.prf.variance"] = flex.pow2(data.intensities.sigmas())
        table.set_flags(flex.bool(table.size(), True), table.flags.integrated_prf)
        correction = flex.double(table.size(), 1.0)
        if data.lp:
            table["lp"] = data.lp
            assert data.lp.all_gt(0.0)
            correction /= data.lp
        if data.qe:
            table["qe"] = data.qe
            correction *= data.qe
        if data.qe or data.lp:
            table["intensity.prf.value"] *= correction
            table["intensity.prf.variance"] *= correction**2

    if data.bg:
        table["background.sum.value"] = data.bg
    if data.bgsig:
        table["background.sum.variance"] = data.bgsig**2
    table["partial_id"] = flex.double(range(0, table.size()))
    table.set_flags(flex.bool(table.size(), True), table.flags.predicted)
    table["id"] = flex.int(table.size(), 0)
    refls_per_scan = OrderedDict({})
    for i, scan in enumerate(scans.values()):
        batch_begin = scan["batch_begin"]
        batch_end = scan["batch_end"]
        sel = (data.batches >= batch_begin) & (data.batches <= batch_end)
        table["id"].set_selected(sel, i)
        refls_per_scan[i + 1] = sel.count(True)
    logger.info(
        "\nNumbers of reflections determined per scan:\n"
        + tabulate(
            [[i, v] for i, v in refls_per_scan.items()],
            ["Scan number", "No. of reflections"],
        )
    )

    # now want to convert phi to z for each sweep
    if params.input.oscillation_step:
        assert params.input.oscillation_step > 0
        z = data.phi / params.input.oscillation_step
    else:
        z = flex.double(table.size(), 0.0)
        angle = data.phi
        for i, expt in enumerate(experiments):
            if expt.scan.get_oscillation()[1] == 0:
                raise ValueError(
                    """Unable to read scan oscillation step.
Please provide an input value for oscillation_step"""
                )
            sel = table["id"] == i
            table.experiment_identifiers()[i] = expt.identifier
            angles_i = angle.select(sel)
            z_i = angles_i / expt.scan.get_oscillation()[1]
            z.set_selected(sel.iselection(), z_i)

    table["xyzobs.px.value"] = flex.vec3_double(data.xpos, data.ypos, z)
    table["xyzcal.px"] = flex.vec3_double(data.xpos, data.ypos, z)
    table["panel"] = flex.size_t(table.size(), 0)
    table["s1"] = flex.vec3_double(table.size(), (0, 0, 0))

    for i, expt in enumerate(experiments):
        fixed_rotation = np.array(expt.goniometer.get_fixed_rotation()).reshape(3, 3)
        setting_rotation = np.array(expt.goniometer.get_setting_rotation()).reshape(
            3, 3
        )
        rotation_axis = np.array(expt.goniometer.get_rotation_axis_datum())
        sel = table["id"] == i
        phi = flumpy.to_numpy(
            expt.scan.get_angle_from_array_index(
                table["xyzobs.px.value"].select(sel).parts()[2], deg=False
            )
        )
        rotation_matrix = Rotation.from_rotvec(
            phi[:, np.newaxis] * rotation_axis
        ).as_matrix()
        R = setting_rotation @ rotation_matrix @ fixed_rotation
        s0 = expt.beam.get_s0()
        UB = np.array(expt.crystal.get_A()).reshape(3, 3)
        total_transform = R @ UB
        hkl = flumpy.to_numpy(table["miller_index"].as_vec3_double())
        s1 = np.einsum("ijk,ik->ij", total_transform, hkl)
        table["s1"].set_selected(sel, flex.vec3_double([row + s0 for row in s1]))

    return table


def create_experiments_from_mtz_without_image_template(
    unmerged_mtz: iotbx.mtz.object,
    params: libtbx.phil.scope_extract,
) -> ExperimentList:

    miller_arrays = unmerged_mtz.as_miller_arrays(
        merge_equivalents=False, anomalous=True
    )
    intensities = None
    for array in miller_arrays:
        if array.info().labels == params.input.labels.intensity.split(","):
            intensities = array
            continue
    if not intensities:
        raise ValueError(
            f"Unable to extract intensities using labels: {params.input.labels.intensity}"
        )

    unit_cell = intensities.unit_cell()
    space_group = intensities.crystal_symmetry().space_group()
    logger.info(
        f"Space group : {space_group.info()} \nUnit cell parameters: {str(unit_cell)}"
    )

    scans = scan_info_from_batch_headers(unmerged_mtz)

    # Get the wavelength for the Beam object
    wavelengths = set()
    for c in unmerged_mtz.crystals():
        for d in c.datasets():
            if d.n_columns() > 0:
                wavelengths.add(d.wavelength())
    if len(wavelengths) > 1:
        raise ValueError(
            f"""More than one wavelength found in the MTZ file: {",".join(wavelengths)},
only single wavelength MTZ files supported."""
        )  ## FIXME allow multi-wave?
    wavelength = list(wavelengths)[0]
    if wavelength == 0.0:
        # try again
        wavelengths = set()
        for c in unmerged_mtz.crystals():
            for d in c.datasets():
                if d.wavelength() > 0.0:
                    wavelengths.add(d.wavelength())
        if len(wavelengths) != 1:
            raise ValueError("Unable to find a single non-zero wavelength")
        wavelength = list(wavelengths)[0]

    # Get detector info from the headers
    panel_size, panel_distance = detector_info_from_batch_headers(unmerged_mtz)

    # Note, for the detector, any refinement of the origin and axes are lost -
    # this information is not output to the MTZ files. But this shouldn't matter
    # for post-integration work.

    experiments = ExperimentList()
    Bmat = mosflm_B_matrix(unit_cell)
    for i, scan in enumerate(scans.values()):
        # we weren't given the images, so have to try and reconstruct as much
        # as possible from the mtz data.
        n_images = scan["end_image"] - (scan["start_image"] - 1)
        if params.input.oscillation_step:
            assert params.input.oscillation_step > 0
            logger.info(
                "Using user-specified oscillation step for constructing the scan"
            )
            dxscan = Scan(
                image_range=[scan["batch_begin"], scan["batch_end"]],
                oscillation=[scan["angle_begin"], params.input.oscillation_step],
            )
        else:
            dxscan = Scan(
                image_range=[scan["batch_begin"], scan["batch_end"]],
                oscillation=[
                    scan["angle_begin"],
                    (scan["angle_end"] - scan["angle_begin"]) / n_images,
                ],
            )

        beam = Beam(
            direction=(0, 0, 1.0), wavelength=wavelength
        )  # we're enforcing DIALS geometry here
        goniometer = GoniometerFactory.from_phil(params.input)  ###FIXME - can
        # we determine this from the mtz - are we only outputting one gonio axis
        # in export_mtz. Perhaps for single sweep, can take several reflections
        # and the UB matrix to solve the goniometer orientation?
        if not goniometer:
            # do the default goniometer
            goniometer = GoniometerFactory.single_axis()
        Umat = matrix.sqr(
            scan["umat_start"]
        ).transpose()  # Convert from Fortranic ordering
        UB = scan["transform_to_dials"] * Umat * Bmat
        ### FIXME do we need to transform by the fixed rotation?
        F = matrix.sqr(goniometer.get_fixed_rotation())
        crystal = Crystal(A=(F.inverse() * UB), space_group=space_group)
        d = Detector()
        p = d.add_panel()
        p.set_image_size(panel_size)
        origin = p.get_origin()
        p.set_frame((1, 0, 0), (0, -1, 0), (origin[0], origin[1], -1 * panel_distance))
        experiments.append(
            Experiment(
                beam=beam,
                scan=dxscan,
                goniometer=goniometer,
                crystal=crystal,
                detector=d,
            )
        )

    set_initial_valid_image_ranges(experiments)
    return experiments


def update_experiments_using_MTZ(
    experiments: ExperimentList,
    unmerged_mtz: iotbx.mtz.object,
    intensity_labels: List[str] = None,
) -> ExperimentList:
    """Add crystals to the experiments using data in the MTZ batch headers."""
    if not intensity_labels:
        intensity_labels = ["I", "SIGI"]

    miller_arrays = unmerged_mtz.as_miller_arrays(
        merge_equivalents=False, anomalous=True
    )
    intensities = None
    for array in miller_arrays:
        if array.info().labels == intensity_labels:
            intensities = array
            continue
    if not intensities:
        raise ValueError(
            f"Unable to extract intensities using labels: {intensity_labels}"
        )

    unit_cell = intensities.unit_cell()
    space_group = intensities.crystal_symmetry().space_group()
    logger.info(
        f"Space group : {space_group.info()} \nUnit cell parameters: {str(unit_cell)}"
    )

    scans = scan_info_from_batch_headers(unmerged_mtz)

    # check that the number of experiments generated from dxtbx matches the
    #  number of scans intepreted from the MTZ
    if len(experiments) != len(scans.values()):
        raise ValueError(
            f"""Mismatch between number of experiments intepreted from image files ({len(experiments)})
and number of scans intepreted from MTZ file ({len(scans.values())}"""
        )
    ## FIXME - what about large sweeps where there is no data at the edge of the
    ## sweep, resulting in no data for those batches in the mtz?
    for i, (expt, scan) in enumerate(zip(experiments, scans.values())):
        if expt.scan.get_image_range() != (scan["start_image"], scan["end_image"]):
            raise ValueError(
                f"""
Inconsistency between image range in experiments and mtz headers for sweep {i+1}.
Scan image range from experiment: {expt.scan.get_image_range()}
Scan image range from mtz headers: {(scan["start_image"], scan["end_image"])}"""
            )

    Bmat = mosflm_B_matrix(unit_cell)
    for i, scan in enumerate(scans.values()):
        Umat = matrix.sqr(
            scan["umat_start"]
        ).transpose()  # Convert from Fortranic ordering
        UB = scan["transform_to_dials"] * Umat * Bmat
        # In dials, we apply F to U before exporting to MTZ, so if we have
        # this we want to use it to get back to the dials UB
        ### FIXME what about mtz from other programs?
        F = matrix.sqr(experiments[0].goniometer.get_fixed_rotation())
        crystal = Crystal(A=(F.inverse() * UB), space_group=space_group)
        experiments[i].crystal = crystal

    return experiments


def create_experiments_from_mtz(
    unmerged_mtz: iotbx.mtz.object,
    params: libtbx.phil.scope_extract,
) -> ExperimentList:
    ### Ideally, the image template can be provided. This allows reading of
    ### the (unrefined) detector models, goniometer, beam and scan. Then, we
    ### only need to fill in the crystal model from the MTZ data. This will
    ### invariably lose some of the scan varying models from refinement, however
    ### these have minimal direct effect on the data reduction process.
    if params.input.images.template:
        experiments = ExperimentList.from_templates(
            [params.input.images.template],
            image_range=params.input.images.image_range,
        )
        experiments = update_experiments_using_MTZ(
            experiments,
            unmerged_mtz,
            intensity_labels=params.input.labels.intensity.split(","),
        )
    else:
        experiments = create_experiments_from_mtz_without_image_template(
            unmerged_mtz, params
        )
    # now set identifiers
    generate_experiment_identifiers(experiments)
    return experiments


def update_experiments_and_reflection_table(
    experiments: ExperimentList, reflections: flex.reflection_table
):
    """Once both an experiment list and a reflection table are available,
    the detector model can be updated by fitting from reflections, and the
    mm position data can be set in reflections. This allows
    dials.two_theta_refine to run on the imported data even if an image
    template is not provided. The accuracy of the resulting geometry cannot
    be guaranteed."""

    # Set pixel size and panel origins
    for i_expt, expt in enumerate(experiments):

        sel = reflections["id"] == i_expt
        ref = reflections.select(reflections["id"] == i_expt)

        for i_panel, panel in enumerate(expt.detector):
            if panel.get_pixel_size() != (0, 0):
                continue

            dn = matrix.col(panel.get_normal())
            df = matrix.col(panel.get_fast_axis())
            ds = matrix.col(panel.get_slow_axis())

            # Calculate cosine of angles between s1 and detector normal
            us0 = matrix.col(expt.beam.get_unit_s0())
            if dn.dot(us0) < 0:
                dn *= -1
            s1 = ref["s1"].select(ref["panel"] == i_panel)
            us1 = s1 / s1.norms()
            cos_tilt = us1.dot(dn)

            # Calculate the scale factor to project s1 onto the detector plane
            scale = panel.get_distance() / cos_tilt
            s1_projected = us1 * scale

            # Get x, y positions in mm from the point of intersection with detector normal
            x_mm = s1_projected.dot(df)
            y_mm = s1_projected.dot(ds)

            # Get x, y positions in pixels
            x_px, y_px, _ = ref["xyzcal.px"].parts()

            # Perform fit to get the origin and pixel size
            x_fit = linregress(x_px, x_mm)
            y_fit = linregress(y_px, y_mm)

            pixel_size_x, pixel_size_y = round(x_fit.slope, 3), round(y_fit.slope, 3)
            shift_x, shift_y = x_fit.intercept, y_fit.intercept
            if pixel_size_x < 0:
                pixel_size_x *= -1
                shift_x *= -1
            if pixel_size_y < 0:
                pixel_size_y *= -1
                shift_y *= -1
            origin = (
                panel.get_distance() * dn + float(shift_x) * df + float(shift_y) * ds
            )

            panel.set_frame(df, ds, origin)
            panel.set_pixel_size((pixel_size_x, pixel_size_y))

    # Update the mm positions
    if "xyzobs.mm.value" not in reflections:
        reflections["xyzcal.mm"] = flex.vec3_double(len(reflections))

        for i_expt, expt in enumerate(experiments):
            sel = reflections["id"] == i_expt
            ref = reflections.select(reflections["id"] == i_expt)
            x, y, z = ref["xyzcal.px"].parts()
            ref["phi"] = expt.scan.get_angle_from_array_index(z)
            ray_intersection(expt.detector, ref)

            reflections["xyzcal.mm"].set_selected(sel, ref["xyzcal.mm"])

        reflections["xyzobs.mm.value"] = reflections["xyzcal.mm"]

    # Set dummy variances to allow refinement
    if "xyzobs.mm.variance" not in reflections:
        reflections["xyzobs.mm.variance"] = flex.vec3_double(len(reflections))

        for i_expt, expt in enumerate(experiments):
            sel = reflections["id"] == i_expt
            ref = reflections.select(reflections["id"] == i_expt)
            px_size = expt.detector[0].get_pixel_size()
            im_width = expt.scan.get_oscillation()[1]
            var_x = flex.double(len(ref), (px_size[0] / 2.0) ** 2)
            var_y = flex.double(len(ref), (px_size[1] / 2.0) ** 2)
            var_phi = flex.double(len(ref), (im_width / 2.0) ** 2)

        reflections["xyzobs.mm.variance"].set_selected(
            sel, flex.vec3_double(var_x, var_y, var_phi)
        )


@show_mail_handle_errors()
def run(args: List[str] = None, phil: libtbx.phil.scope = phil_scope) -> None:
    """Run the command-line script."""

    usage = "dials.import_mtz integrated.mtz"

    parser = OptionParser(
        usage=usage,
        phil=phil,
        epilog=help_message,
        read_experiments=False,
        read_reflections=False,
        check_format=False,
    )

    params, _, unhandled = parser.parse_args(
        args=args, return_unhandled=True, show_diff_phil=True
    )

    log.config(logfile="dials.import_mtz.log")

    if len(unhandled) == 1:
        mtzfile = unhandled[0]
    else:
        raise ValueError("Only one input MTZ file expected")

    logger.info(f"\nReading data from MTZ file: {mtzfile}\n")
    unmerged_mtz = iotbx.mtz.object(file_name=mtzfile)

    if len(params.input.labels.intensity.split(",")) != 2:
        sys.exit(
            str(
                ValueError(
                    "Intensity labels given must be two comma-separated values e.g 'I,SIGI' "
                )
            )
        )

    try:
        experiments = create_experiments_from_mtz(unmerged_mtz, params)
        table = create_reflection_table_from_mtz(unmerged_mtz, params, experiments)
    except ValueError as e:
        sys.exit(str(e))

    update_experiments_and_reflection_table(experiments, table)

    logger.info(f"Saving extracted data to {params.output.reflections}")
    table.as_file(params.output.reflections)

    logger.info(f"Saving extracted metadata to {params.output.experiments}")
    experiments.as_file(params.output.experiments)


if __name__ == "__main__":
    run()

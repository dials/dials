"""
dials.import_mtz, a tool for creating DIALS working data formats (reflection
tables and experiments) from an MTZ file.
"""

from __future__ import absolute_import, division, print_function

import logging
import math
import sys
from collections import OrderedDict

from dxtbx.model import (
    Beam,
    Crystal,
    Detector,
    Experiment,
    ExperimentList,
    Goniometer,
    Scan,
)
from iotbx import mtz
from libtbx.phil import parse
from scitbx import matrix

from dials.algorithms.refinement.refinement_helpers import set_obs_s1
from dials.array_family import flex
from dials.util import log, show_mail_handle_errors, tabulate
from dials.util.options import OptionParser

logger = logging.getLogger("dials.command_line.import_mtz")

help_message = """

This program makes DIALS datafiles from an MTZ input.
"""

phil_scope = parse(
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
        pixel_size = 0.172
            .type = float(value_min=0.0)
    }

  }

  output {

    reflections = imported_mtz.refl
      .type = str
      .help = "Name of the output reflections file"

    experiments = imported_mtz.expt
      .type = str
      .help = "Name of the output experiments file"

  }
"""
)


def detector_info_from_batch_headers(unmerged_mtz):
    b0 = unmerged_mtz.batches()[0]
    panel_size = int(b0.detlm()[1]), int(b0.detlm()[3])
    panel_distance = b0.dx()[0]
    return panel_size, panel_distance


def mosflm_B_matrix(unit_cell):
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


def scan_info_from_batch_headers(unmerged_mtz, verbose=True):
    batches = unmerged_mtz.batches()

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
            }
        }
    )

    scan_no = 1
    phi_end = batches[0].phiend()
    last_batch = batches[0].num()

    for b in batches[1:]:
        if abs(b.phistt() - phi_end) > 0.0001:
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
                "umat_start": b.umat(),
                "s0n": b.so(),
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


def create_reflection_table_from_mtz(
    unmerged_mtz, params, experiments, batch_offset=None
):
    """Extract reflection table data"""

    miller_arrays = unmerged_mtz.as_miller_arrays(
        merge_equivalents=False, anomalous=True
    )

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

    # Select the desired columns
    intensities = None
    profile_intensities = None
    batches = None
    xpos = None
    ypos = None
    phi = None
    scale_used = None
    # optional columns?
    partiality = None
    bg = None
    bgsig = None
    lp = None
    qe = None
    columns_found = []
    for array in miller_arrays:
        labels = array.info().labels
        columns_found.extend(labels)
        if labels == [params.input.labels.xpos]:
            xpos = array
        elif labels == [params.input.labels.ypos]:
            ypos = array
        elif labels == [params.input.labels.phi]:
            phi = array  # this is xyzcal.px
        elif labels == [params.input.labels.partiality]:
            partiality = array
            # want to transform from rotation angle to image number
        elif labels == params.input.labels.intensity.split(","):
            intensities = array
        elif labels == ["IPR", "SIGIPR"]:
            profile_intensities = array
        elif labels == ["BATCH"]:
            batches = array
        elif labels == ["BG"]:
            bg = array
        elif labels == ["SIGBG"]:
            bgsig = array
        elif labels == ["LP"]:
            lp = array
        elif labels == ["QE"]:
            qe = array
        elif labels == ["SCALEUSED"]:
            scale_used = array
    logger.info(f"Data arrays found in mtz file: {', '.join(columns_found)}")
    if not intensities:
        raise ValueError(
            f"Intensities not found using labels {params.input.labels.intensity}"
            + "\nIntensity column labels can be specified with labels.intensity='Iname,SIGIname'"
        )
    if not batches:
        raise ValueError("Batch values not found")
    if not phi:
        raise ValueError(
            f"""Phi/Rot values not found using labels {params.input.labels.phi}
Phi/Rot column label can be specified with labels.phi='PHIname'"""
        )
    if batches.data().size() != intensities.data().size():
        raise ValueError("Batch and intensity array sizes do not match")
    batches = batches.data()
    # for some reason the indices aren't actually the indices, so do this:
    indices = unmerged_mtz.extract_original_index_miller_indices()

    scans = scan_info_from_batch_headers(unmerged_mtz, verbose=False)

    # The reflection data
    table = flex.reflection_table()
    table["miller_index"] = indices
    table["d"] = intensities.d_spacings().data()
    if partiality:
        table["partiality"] = partiality.data()
    else:
        table["partiality"] = flex.double(table.size(), 1.0)

    # extract the intensities. In dials, we work on raw intensities.
    # In mtz files, the intensities already have corrections applied. i.e
    # if we have IPR and LP, then intensity.prf.value = IPR/LP

    # First, if the SCALEUSED column is there, indicates scaled data
    # In that case, I,SIGI are scaled intensities
    if scale_used:
        logger.info(
            """The intensity data with labels I,SIGI are interpeted as being
scaled intensity data. Any further scaling correction will be made in
addition to the existing correction (SCALEUSED)."""
        )
        # Assertion is that we want to use the scaled intensities.
        # Set the scaled intensities as 'intensity.prf.value'
        table["intensity.prf.value"] = intensities.data()
        table["intensity.prf.variance"] = flex.pow2(intensities.sigmas())
        table["scaledused"] = scale_used.data()
        # now these have had corrections like lp applied already. So won't
        # be wanting to apply these again, therefore don't set these columns.
        table.set_flags(flex.bool(table.size(), True), table.flags.integrated_prf)

    # FIXME check these assumptions for aimless, xds, etc data
    # now see whether two or one sets of intensities were found
    elif profile_intensities and intensities:
        logger.info(
            f"""The intensity data with labels I,SIGI are interpeted as being
summation integrated intensities.
The data with labels IPR,SIGIPR are interpreted as being
profile-fitted integrated intensities."""
        )
        table["intensity.sum.value"] = intensities.data()
        table["intensity.sum.variance"] = flex.pow2(intensities.sigmas())
        table["intensity.prf.value"] = profile_intensities.data()
        table["intensity.prf.variance"] = flex.pow2(profile_intensities.sigmas())
        table.set_flags(flex.bool(table.size(), True), table.flags.integrated_sum)
        table.set_flags(flex.bool(table.size(), True), table.flags.integrated_prf)
        correction = flex.double(table.size(), 1.0)
        if lp:
            table["lp"] = lp.data()
            assert lp.data().all_gt(0.0)
            correction /= lp.data()
        if qe:
            table["qe"] = qe.data()
            correction *= qe.data()
        if qe or lp:
            table["intensity.sum.value"] *= correction
            table["intensity.sum.variance"] *= correction ** 2
            table["intensity.prf.value"] *= correction
            table["intensity.prf.variance"] *= correction ** 2
        if partiality:
            table["intensity.sum.value"] *= partiality.data()
            table["intensity.sum.variance"] *= partiality.data() ** 2

    else:  # only one intensity set, assume profile intensities.
        logger.info(
            """The data with labels I,SIGI are interpreted as being
profile-fitted integrated intensities."""
        )
        table["intensity.prf.value"] = intensities.data()
        table["intensity.prf.variance"] = flex.pow2(intensities.sigmas())
        table.set_flags(flex.bool(table.size(), True), table.flags.integrated_prf)
        correction = flex.double(table.size(), 1.0)
        if lp:
            table["lp"] = lp.data()
            assert lp.data().all_gt(0.0)
            correction /= lp.data()
        if qe:
            table["qe"] = qe.data()
            correction *= qe.data()
        if qe or lp:
            table["intensity.prf.value"] *= correction
            table["intensity.prf.variance"] *= correction ** 2

    if bg:
        table["background.sum.value"] = bg.data()
    if bgsig:
        table["background.sum.variance"] = bgsig.data() ** 2
    table["partial_id"] = flex.double(range(0, table.size()))

    table["id"] = flex.int(table.size(), 0)
    refls_per_scan = OrderedDict({})
    for i, scan in enumerate(scans.values()):
        batch_begin = scan["batch_begin"]
        batch_end = scan["batch_end"]
        sel = (batches >= batch_begin) & (batches <= batch_end)
        table["id"].set_selected(sel, i)
        refls_per_scan[i + 1] = sel.count(True)
    logger.info(
        "Numbers of reflections determined per scan:\n"
        + tabulate(
            [[i, v] for i, v in refls_per_scan.items()],
            ["Scan number", "No. of reflections"],
        )
    )

    # now want to convert phi to z for each sweep
    if params.input.oscillation_step:
        assert params.input.oscillation_step > 0
        z = phi.data() / params.input.oscillation_step
    else:
        z = flex.double(table.size(), 0.0)
        angle = phi.data()
        for i, expt in enumerate(experiments):
            if expt.scan.get_oscillation()[1] == 0:
                raise ValueError(
                    """Unable to read scan oscillation step.
Please provide an input value for oscillation_step"""
                )
            sel = table["id"] == 0
            angles_i = angle.select(sel)
            z_i = angles_i / expt.scan.get_oscillation()[1]
            z.set_selected(sel.iselection(), z_i)

    table["xyzobs.px.value"] = flex.vec3_double(xpos.data(), ypos.data(), z)
    table["xyzcal.px"] = flex.vec3_double(xpos.data(), ypos.data(), z)
    table["panel"] = flex.size_t(table.size(), 0)

    table.centroid_px_to_mm(experiments)
    table["s1"] = flex.vec3_double(table.size(), (0, 0, 0))
    set_obs_s1(table, experiments)

    return table


def create_experiments_from_MTZ(unmerged_mtz, params):

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
        Umat = matrix.sqr(scan["umat_start"])
        ### FIXME do we need to transform by the fixed rotation?
        crystal = Crystal(A=Umat.transpose() * Bmat, space_group=space_group)

        # we weren't given the images, so have to try and reconstruct as much
        # as possible from the mtz data.
        n_images = scan["end_image"] - (scan["start_image"] - 1)
        if params.input.oscillation_step:
            assert params.input.oscillation_step > 0
            logger.info(
                "Using user-specified oscillation step for constructing the scan"
            )
            dxscan = Scan(
                image_range=[scan["start_image"], scan["end_image"]],
                oscillation=[0.0, params.input.oscillation_step],
            )
        else:
            dxscan = Scan(
                image_range=[scan["start_image"], scan["end_image"]],
                oscillation=[
                    0.0,
                    (scan["angle_end"] - scan["angle_begin"]) / n_images,
                ],
            )

        beam = Beam(direction=tuple(-1.0 * scan["s0n"]), wavelength=wavelength)
        goniometer = Goniometer(
            (1.0, 0.0, 0.0)
        )  ###FIXME - we are only outputting one gonio axis in export_mtz
        d = Detector()
        p = d.add_panel()
        p.set_image_size(panel_size)
        p.set_pixel_size(
            (params.input.detector.pixel_size, params.input.detector.pixel_size)
        )
        experiments.append(
            Experiment(
                beam=beam,
                scan=dxscan,
                goniometer=goniometer,  ##FIXME
                crystal=crystal,  ## FIXME - consistent UB settings
                detector=d,  ##FIXME - do we need to allow setting more detector parameters?
            )
        )

    return experiments


def update_experiments_using_MTZ(experiments, unmerged_mtz, intensity_labels=None):
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
        Umat = matrix.sqr(scan["umat_start"])
        # In dials, we apply F to U before exporting to MTZ, so if we have
        # this we want to use it to get back to the dials UB
        ### FIXME what about mtz from other programs?
        F = matrix.sqr(experiments[0].goniometer.get_fixed_rotation())
        crystal = Crystal(
            A=(F.inverse() * Umat.transpose() * Bmat), space_group=space_group
        )
        experiments[i].crystal = crystal

    return experiments


@show_mail_handle_errors()
def run(args=None, phil=phil_scope):
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
    unmerged_mtz = mtz.object(file_name=mtzfile)

    if len(params.input.labels.intensity.split(",")) != 2:
        sys.exit(
            str(
                ValueError(
                    "Intensity labels given must be two comma-separated values e.g 'I,SIGI' "
                )
            )
        )

    ### Ideally, the image template can be provided. This allows reading of
    ### the (unrefined) detector models, goniometer, beam and scan. Then, we
    ### only need to fill in the crystal model from the MTZ data. This will
    ### invariably lose some of the scan varying models from refinement, however
    ### these have minimal direct effect on the data reduction process.
    if params.input.images.template:
        try:
            experiments = ExperimentList.from_templates(
                [params.input.images.template],
                image_range=params.input.images.image_range,
            )
            experiments = update_experiments_using_MTZ(
                experiments,
                unmerged_mtz,
                intensity_labels=params.input.labels.intensity.split(","),
            )
        except ValueError as e:
            sys.exit(str(e))
    else:
        experiments = create_experiments_from_MTZ(unmerged_mtz, params)

    try:
        table = create_reflection_table_from_mtz(unmerged_mtz, params, experiments)
    except ValueError as e:
        sys.exit(str(e))

    logger.info(f"Saving extracted data to {params.output.reflections}")
    table.as_file(params.output.reflections)

    logger.info(f"Saving extracted metadata to {params.output.experiments}")
    experiments.as_file(params.output.experiments)


if __name__ == "__main__":
    run()

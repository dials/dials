from __future__ import absolute_import, division, print_function

import logging
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
from dxtbx.model.experiment_list import (
    ExperimentListFactory,
    ExperimentListTemplateImporter,
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
            .multiple = True

        directory = None
            .type = str
            .help = "A directory with images"
            .multiple = True

        image = None
            .type = str
            .help = "An image for the dataset"
            .multiple = True
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


def scan_info_from_batch_headers(unmerged_mtz):
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

    return scans


def read_mtzfile(filename, params, experiments=None, batch_offset=None):
    """
    Read the mtz file
    """
    unmerged_mtz = mtz.object(file_name=filename)
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
    # optional columns?
    partiality = None
    bg = None
    bgsig = None
    lp = None
    qe = None
    for array in miller_arrays:
        if array.info().labels == [params.input.labels.xpos]:
            xpos = array
        elif array.info().labels == [params.input.labels.ypos]:
            ypos = array
        elif array.info().labels == [params.input.labels.phi]:
            phi = array  # this is xyzcal.px
        elif array.info().labels == [params.input.labels.partiality]:
            partiality = array
            # want to transform from rotation angle to image number
        elif array.info().labels == ["I", "SIGI"]:
            intensities = array
        elif array.info().labels == ["IPR", "SIGIPR"]:
            profile_intensities = array
        elif array.info().labels == ["BATCH"]:
            batches = array
        elif array.info().labels == ["BG"]:
            bg = array
        elif array.info().labels == ["SIGBG"]:
            bgsig = array
        elif array.info().labels == ["LP"]:
            lp = array
        elif array.info().labels == ["QE"]:
            qe = array
    if not intensities:
        raise KeyError("Intensities not found in mtz file, expected labels I, SIGI")
    if not batches:
        raise KeyError("Batch values not found")
    if not phi:
        raise KeyError(
            "Phi/Rot values not found. Column label can be specified with labels.phi="
        )
    if batches.data().size() != intensities.data().size():
        raise ValueError("Batch and intensity array sizes do not match")
    batches = batches.data()

    # now sort out what the intensity columns are.
    # if

    logger.info("Reading batch headers to determine scan properties")
    scans = scan_info_from_batch_headers(unmerged_mtz)
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

    filled_experiments = None
    if not experiments:
        experiments = ExperimentList()
        panel_size, panel_distance = detector_info_from_batch_headers(unmerged_mtz)
    else:
        # check that the number of experiments generated from dxtbx matches the
        #  number of scans intepreted from the MTZ
        nscans = len(scans.values())
        if len(experiments) != nscans:
            raise ValueError(
                f"""
Mismatch between number of experiments intepreted from image files ({len(experiments)})
and number of scans intepreted from MTZ file ({nscans}"""
            )
        filled_experiments = experiments

    # Get the unit cell and space group
    unit_cell = intensities.unit_cell()
    space_group = intensities.crystal_symmetry().space_group()
    logger.info(f"Space group : {space_group.info()} \nUnit cell : {str(unit_cell)}")
    Bmat = matrix.sqr(unit_cell.fractionalization_matrix()).transpose()

    for i, scan in enumerate(scans.values()):
        Umat = matrix.sqr(scan["umat_start"])
        crystal = Crystal(A=Umat * Bmat, space_group=space_group)
        if filled_experiments:
            experiments[i].crystal = crystal
        else:
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
            beam = Beam(s0=tuple(list(scan["s0n"])))
            goniometer = Goniometer((1.0, 0.0, 0.0))  ###FIXME
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
                    goniometer=goniometer,
                    crystal=crystal,
                    detector=d,
                )
            )
    if not filled_experiments:
        filled_experiments = experiments

    # for some reason the indices aren't actually the indices, so do this:
    indices = mtz.object(file_name=filename).extract_original_index_miller_indices()
    # i_obs = i_obs.customized_copy(indices=indices, info=i_obs.info())

    # The reflection data
    table = flex.reflection_table()
    table["miller_index"] = indices  # intensities.indices()
    table["d"] = intensities.d_spacings().data()

    # now see whether two or one sets of intensities were found
    if profile_intensities:
        table["intensity.sum.value"] = intensities.data()
        table["intensity.sum.variance"] = flex.pow2(intensities.sigmas())
        table["intensity.prf.value"] = profile_intensities.data()
        table["intensity.prf.variance"] = flex.pow2(profile_intensities.sigmas())
        table.set_flags(flex.bool(table.size(), True), table.flags.integrated_sum)
        table.set_flags(flex.bool(table.size(), True), table.flags.integrated_prf)
    else:  # only one intensity set, assume profile intensities.
        table["intensity.prf.value"] = intensities.data()
        table["intensity.prf.variance"] = flex.pow2(intensities.sigmas())
        table.set_flags(flex.bool(table.size(), True), table.flags.integrated_prf)
    if partiality:
        table["partiality_applied"] = partiality.data()
    if bg:
        table["background.sum.value"] = bg.data()
    if bgsig:
        table["background.sum.variance"] = bgsig.data() ** 2
    if lp:
        table["lp_applied"] = lp.data()
    if qe:
        table["qe_applied"] = qe.data()
    table["partial_id"] = flex.double(range(0, table.size()))

    # need s1 ##FIXME
    # table["s1"] = flex.vec3_double(table.size(), (1, 1, 1))

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
                    """
    Unable to read scan oscillation step.
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

    return table, experiments


def load_from_images(params):
    experiments = None
    if params.input.images.directory:

        experiments = ExperimentListFactory.from_filenames(
            params.input.images.directory
        )
    elif params.input.images.template:
        experiments = ExperimentListTemplateImporter(
            params.input.images.template
        ).experiments
    elif params.input.images.image:
        from dials.util.options import Importer, flatten_experiments

        experiments = flatten_experiments(
            Importer(
                [params.input.images.image],
                read_experiments=False,
                read_reflections=False,
                read_experiments_from_images=True,
            ).experiments
        )
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

    log.config(logfile="dials.command_line.import_mtz")

    if len(unhandled) == 1:
        mtzfile = unhandled[0]
    else:
        raise ValueError("Only one input MTZ file expected")

    try:
        experiments = load_from_images(params)
    except ValueError as e:
        sys.exit(str(e))

    try:
        table, expts = read_mtzfile(mtzfile, params, experiments=experiments)
    except ValueError as e:
        sys.exit(str(e))

    logger.info(f"Saving extracted data to {params.output.reflections}")
    table.as_file(params.output.reflections)

    logger.info(f"Saving extracted metadata to {params.output.experiments}")
    expts.as_file(params.output.experiments)


if __name__ == "__main__":
    run()

from __future__ import absolute_import, division, print_function

import logging
from collections import OrderedDict

from iotbx import mtz
from libtbx.phil import parse

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

    mtzfile = None
      .type = str
      .help = "We can also import an MTZ file"

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

  }
"""
)


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
            }

        phi_end = b.phiend()
        last_batch = b.num()

    scans[scan_no]["angle_end"] = phi_end
    scans[scan_no]["batch_end"] = last_batch
    scans[scan_no]["end_image"] = last_batch - scans[scan_no]["batch_begin"] + 1

    return scans


def read_mtzfile(filename, params, batch_offset=None):
    """
    Read the mtz file
    """
    miller_arrays = mtz.object(file_name=filename).as_miller_arrays(
        merge_equivalents=False, anomalous=True
    )

    # Select the desired columns
    intensities = None
    batches = None
    xpos = None
    ypos = None
    phi = None
    partiality = None
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
        elif array.info().labels == ["BATCH"]:
            batches = array
    if not intensities:
        raise KeyError("Intensities not found in mtz file, expected labels I, SIGI")
    if not batches:
        raise KeyError("Batch values not found")
    if batches.data().size() != intensities.data().size():
        raise ValueError("Batch and intensity array sizes do not match")
    batches = batches.data()

    logger.info("Reading batch headers to determine scan properties")
    scans = scan_info_from_batch_headers(mtz.object(file_name=filename))
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

    # Get the unit cell and space group
    unit_cell = intensities.unit_cell()
    space_group = intensities.crystal_symmetry().space_group()
    logger.info(f"Space group : {space_group.info()} \nUnit cell : {str(unit_cell)}")

    # for some reason the indices aren't actually the indices, so do this:
    indices = mtz.object(file_name=filename).extract_original_index_miller_indices()
    # i_obs = i_obs.customized_copy(indices=indices, info=i_obs.info())

    # The reflection data
    table = flex.reflection_table()
    table["miller_index"] = indices  # intensities.indices()
    table["d"] = intensities.d_spacings().data()
    table["intensity"] = intensities.data()
    table["variance"] = flex.pow2(intensities.sigmas())
    table["partiality"] = partiality.data()
    table["partial_id"] = flex.double(range(0, table.size()))

    # need s1
    table["s1"] = flex.vec3_double(table.size(), (1, 1, 1))

    table["id"] = flex.int(table.size(), 0)
    refls_per_scan = OrderedDict({})
    for i, scan in enumerate(scans.values()):
        batch_begin = scan["batch_begin"]
        batch_end = scan["batch_end"]
        sel = (batches >= batch_begin) & (batches <= batch_end)
        table["id"].set_selected(sel, i)
        refls_per_scan[i + 1] = sel.count(True)
    logger.info("Numbers of reflections determined per scan:")
    logger.info(
        tabulate(
            [[i, v] for i, v in refls_per_scan.items()],
            ["Scan number", "No. of reflections"],
        )
    )

    # now want to convert phi to z for each sweep
    if not phi:
        z = flex.double(table.size(), 0)
        for id_ in set(table["id"]):
            sel = table["id"] == id_
            phi_i = phi.select(sel)
            if params.input.oscillation_step:
                assert params.input.oscillation_step > 0
                z_i = (phi_i.data() - min(phi_i.data())) / params.input.oscillation_step
            else:
                # try to approximate using batch
                phirange = int(flex.max(phi_i.data())) - int(flex.min(phi_i.data()))
                # add 1 because a 360 image sweep will have batches 1 to 360
                batches_i = batches.select(sel)
                batch_range = max(batches_i) - min(batches_i) + 1
                z_i = (phi_i.data() - flex.min(phi_i.data())) * batch_range / phirange
                print(batch_range, phirange, max(z_i), min(z_i))
            z.set_selected(sel.iselection(), z_i)
    else:
        # found a phi column
        if params.input.oscillation_step:
            assert params.input.oscillation_step > 0
            z = phi.data() / params.input.oscillation_step
        else:
            z = phi.data()  # use scans info to determine step.

    table["xyzobs.px.value"] = flex.vec3_double(xpos.data(), ypos.data(), z)

    return table, unit_cell, space_group


@show_mail_handle_errors()
def run(args=None, phil=phil_scope):
    """Run the command-line script."""

    usage = "dials.import_mtz input.mtzfile=example.mtz"

    parser = OptionParser(
        usage=usage,
        phil=phil,
        epilog=help_message,
        read_experiments=False,
        read_reflections=False,
        check_format=False,
    )

    params, _ = parser.parse_args(args=args, show_diff_phil=False)

    log.config(logfile="dials.command_line.import_mtz")

    table, unit_cell, space_group = read_mtzfile(params.input.mtzfile, params)

    output = "test.refl"
    logger.info(f"Saving extracted data to {output}")
    table.as_file(output)


if __name__ == "__main__":
    run()

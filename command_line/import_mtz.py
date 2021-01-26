from __future__ import absolute_import, division, print_function

import logging

from iotbx import mtz
from libtbx.phil import parse

from dials.array_family import flex
from dials.util import log, show_mail_handle_errors
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


def read_mtzfile(filename, params, batch_offset=None):
    """
    Read the mtz file
    """
    miller_arrays = mtz.object(file_name=filename).as_miller_arrays(
        merge_equivalents=False, anomalous=True
    )

    # print(dir(mtz.object(file_name=filename)))
    batch_objs = mtz.object(file_name=filename).batches()

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
            phi = array
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

    for i, batch in enumerate(batch_objs):
        if i < 5:
            print(list(batch.cell()))
            print(dir(batch))

    # Get the unit cell and space group
    unit_cell = intensities.unit_cell()
    space_group = intensities.crystal_symmetry().space_group()

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

    # Create unit cell list
    zeroed_batches = batches.data() - flex.min(batches.data())
    dataset = flex.int(table.size(), 0)
    sorted_batches = flex.sorted(zeroed_batches)
    sel_perm = flex.sort_permutation(zeroed_batches)

    batch_starts = [sorted_batches[0]]
    if not batch_offset:
        previous = 0
        potential_batch_offsets = flex.double()
        for i, b in enumerate(sorted_batches):
            if b - previous > 1:
                potential_batch_offsets.append(b - previous)
                batch_starts.append(b)
            previous = b
        potential = flex.sorted(potential_batch_offsets)
        # potential is a list of low numbers (where images may not have any spots)
        # and larger numbers between batches.
        print(list(potential))
        if len(potential) == 1:
            batch_offset = potential[0]
            logger.info(
                """
Using a batch offset of %s to split datasets.
Batch offset can be specified with mtz.batch_offset=
""",
                batch_offset,
            )
        elif len(potential) > 1:
            diffs = flex.double(
                [potential[i + 1] - p for i, p in enumerate(potential[:-1])]
            )
            i = flex.sort_permutation(diffs)[-1]
            batch_offset = int(potential[i + 1] - (0.2 * diffs[i]))
            logger.info(
                """
Using an approximate batch offset of %s to split datasets.
Batch offset can be specified with mtz.batch_offset=
""",
                batch_offset,
            )
        else:
            batch_offset = 1

    previous = 0
    dataset_no = 0
    for i, b in enumerate(sorted_batches):
        if b - previous > batch_offset:
            dataset_no += 1
        dataset[i] = dataset_no
        previous = b

    table["id"] = flex.int(table.size(), 0)
    table["id"].set_selected(sel_perm, dataset)
    for id_ in set(table["id"]):
        table.experiment_identifiers()[id_] = str(id_)

    # now want to convert phi to z for each sweep
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
            batches_i = batches.select(sel).data()
            batch_range = max(batches_i) - min(batches_i) + 1
            z_i = (phi_i.data() - flex.min(phi_i.data())) * batch_range / phirange
            print(batch_range, phirange, max(z_i), min(z_i))
        z.set_selected(sel.iselection(), z_i)

    table["xyzobs.px.value"] = flex.vec3_double(xpos.data(), ypos.data(), z)
    print(set(table["id"]))
    print(list(batch_starts))
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

    table.as_file("test.refl")

    print(table)
    print(unit_cell)
    print(space_group)


if __name__ == "__main__":
    run()

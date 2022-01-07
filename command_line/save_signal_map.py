# LIBTBX_SET_DISPATCHER_NAME dev.dials.save_signal_map

import h5py
import hdf5plugin
import numpy

import iotbx.phil
from libtbx.phil import parse

import dials.util.masking
from dials.algorithms.spot_finding.factory import SpotFinderFactory
from dials.algorithms.spot_finding.factory import phil_scope as spot_phil
from dials.util import show_mail_handle_errors
from dials.util.options import OptionParser, flatten_experiments

help_message = """

Examples::

  dev.dials.save_signal_map data_1.nxs

"""

phil_scope = iotbx.phil.parse(
    """
masking {
  include scope dials.util.masking.phil_scope
}

output {
  filename = "signal_map.h5"
    .type = path
    .help = "Filename to save signal map"
}

""",
    process_includes=True,
)


@show_mail_handle_errors()
def run(args=None):
    usage = "dev.dials.save_signal_map [options] data_1.nxs"

    parser = OptionParser(
        usage=usage,
        phil=phil_scope,
        read_experiments=True,
        read_experiments_from_images=True,
        epilog=help_message,
    )

    params, options = parser.parse_args(args, show_diff_phil=True)

    # Ensure we have either a data block or an experiment list
    experiments = flatten_experiments(params.input.experiments)
    imagesets = experiments.imagesets()

    assert len(imagesets) == 1

    imageset = imagesets[0]

    first, last = imageset.get_scan().get_array_range()

    mask_params = phil_scope.fetch(parse("")).extract().masking

    detector = imageset.get_detector()

    # Only working with single panel detector for now
    assert len(detector) == 1
    panel = detector[0]
    all_mask = dials.util.masking.generate_mask(imageset, mask_params)[0]

    spot_params = spot_phil.fetch(source=parse("")).extract()
    threshold_function = SpotFinderFactory.configure_threshold(spot_params)

    with h5py.File(params.output.filename, "x") as fout:
        frames = last - first + 1
        fast, slow = panel.get_image_size()

        dset = fout.create_dataset(
            "data",
            (frames, slow, fast),
            chunks=(1, slow, fast),
            **hdf5plugin.Bitshuffle(),
            dtype=numpy.uint8,
        )

        for indx in range(first, last):
            print(indx)
            imageset_mask = imageset.get_mask(indx)[0]
            mask = imageset_mask & all_mask
            data = imageset.get_raw_data(indx)[0].as_double()
            signal_pixels = threshold_function.compute_threshold(data, mask)
            signal = signal_pixels.as_numpy_array().astype(numpy.uint8)
            dset[indx, :, :] = signal


if __name__ == "__main__":
    run()

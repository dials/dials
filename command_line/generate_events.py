# LIBTBX_SET_DISPATCHER_NAME dev.dials.generate_events
# DIALS_ENABLE_COMMAND_LINE_COMPLETION

from dials.util.options import OptionParser
from dials.util.options import flatten_experiments
from dials.util import show_mail_on_error
from dials.array_family import flex
import iotbx.phil

import dxtbx  # noqa: F401; import dependency to find HDF5 library
import h5py
import sys

"""

"""

# Create phil parameters
phil_scope = iotbx.phil.parse(
    """
        input {
          mask = None
            .type = str
            .help = "Bad pixels mask"
          image_range = None
            .type = ints(value_min = 0, size = 2)
            .help = "Image range e.g. 1,1000"
        }

        output {
          events = None
            .type = str
            .help = "The output HDF5 file name"
        }
        """
)


def run():
    # Parse command line arguments
    parser = OptionParser(
        phil=phil_scope, read_experiments=True, read_experiments_from_images=True
    )
    params, options = parser.parse_args(show_diff_phil=True)

    if params.output.events:
        fout = h5py.File(params.output.events, "x")
    else:
        print("No output file created\n")

    experiments = flatten_experiments(params.input.experiments)
    if len(experiments) != 1:
        parser.print_help()
        sys.exit("Please pass images\n")
        return

    imagesets = experiments.imagesets()
    if len(imagesets) != 1:
        sys.exit("Please pass only one set of images\n")
        return

    imageset = imagesets[0]

    if params.input.mask:
        with h5py.File(params.input.mask, "r") as fh:
            bad_pixels = flex.bool(fh["mask"][()])

    if params.input.image_range:
        image_range = params.input.image_range
    else:
        image_range = 0, len(imageset)

    for j in range(*image_range):
        img = imageset.get_raw_data(j)[0]
        msk = imageset.get_mask(j)

        img.set_selected(~msk, -1)
        if params.input.mask:
            img.set_selected(~bad_pixels, -2)

    fout.close()


if __name__ == "__main__":
    with show_mail_on_error():
        run()

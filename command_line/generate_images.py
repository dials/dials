# LIBTBX_SET_DISPATCHER_NAME dev.dials.generate_events
# DIALS_ENABLE_COMMAND_LINE_COMPLETION

"""
Add doc
"""

import dxtbx  # noqa: F401; import dependency to find HDF5 library
import h5py

import os
import argparse

from libtbx import phil
from dials.util import show_mail_on_error

from dials.algorithms.event_mode.events_to_images import events_to_images


class CheckFileExt(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        fin, fin_ext = os.path.splitext(values[0])
        if fin_ext != ".nxs":
            parser.error("Please pass the NeXus-like .nxs file containing events data.")

        condition = any("image_file" in v for v in values)
        if len(values) > 1 and condition is True:
            i = ["image_file" in v for v in values].index(True)
            fout, fout_ext = os.path.splitext(values[i])
            if fout_ext != ".h5":
                print(
                    f"You specified an invalid file extension {fout_ext} for the output "
                    "image file.\n"
                    f"The output images will be saved to {fout}.h5 instead."
                )
                values[i] = f"{fout}.h5"
        setattr(namespace, self.dest, values)


# Create phil parameters
phil_scope = phil.parse(
    """
        input {
          event_data = None
            .type = path
            .multiple = True
            .help = "Input event mode data from timepix detector, Nexus file."
          exposure_time = 0.1
            .type = float
            .help = "Exposure time of each frame."
        }

        output {
          image_file = None
            .type = path
            .help = "The output HDF5 file."
        }
    """
)


def runme():
    # Parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--debug", action="store_const", const=True)
    parser.add_argument("phil", nargs="+", action=CheckFileExt)

    args = parser.parse_args()
    print(args, args.phil)

    cl = phil_scope.command_line_argument_interpreter()
    working_phil = phil_scope.fetch(cl.process_and_fetch(args.phil))
    params = working_phil.extract()

    # If output file directory doesn't yet exist, create it.
    wdir = os.path.dirname(params.output.image_file)
    if not os.path.exists(wdir):
        os.mkdir(wdir)

    with h5py.File(params.input.event_data[0], "r") as f, h5py.File(
        params.output.image_file, "x"
    ) as g:
        events_to_images(f, params.input.exposure_time, g)


if __name__ == "__main__":
    with show_mail_on_error():
        runme()

# LIBTBX_SET_DISPATCHER_NAME dev.dials.generate_events
# DIALS_ENABLE_COMMAND_LINE_COMPLETION

"""
Add doc
"""

import os
import argparse

from libtbx import phil
from dials.util import show_mail_on_error

from dials.algorithms.event_mode.image_to_events import images_to_events


class CheckFileExt(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        if any("mask" in v for v in values):
            i = ["mask" in v for v in values].index(True)
            msk, msk_ext = os.path.splitext(values[i])
            if msk_ext != ".h5":
                print(
                    f"You specified a mask file with an invalid extension {msk_ext}.\n"
                    f"Please convert file to HDF5.\n"
                    f"The program will run anyway without applying the mask.\n"
                )
                values[i] = None
        if any("events" in v for v in values):
            i = ["events" in v for v in values].index(True)
            fout, fout_ext = os.path.splitext(values[i])
            if fout_ext != ".h5":
                print(
                    f"You specified an invalid file extension {fout_ext} for the output "
                    "event file.\n"
                    f"The output events will be saved to {fout}.h5 instead."
                )
                values[i] = f"{fout}.h5"
        else:
            parser.error("Please specify a .h5 output file.")
        setattr(namespace, self.dest, values)


# Create phil parameters
phil_scope = phil.parse(
    """
        input {
          experiments = None
            .type = path
            .multiple = True
            .help = "Input image data from experiment."

          image_range = None
            .type = ints(value_min = 0, size = 2)
            .help = "Define image range tto turn into events e.g. 1,1000"

          exposure_time = 1000.
            .type = float
            .help = "Exposure time of each frame."

          mask = None
            .type = path
            .multiple = True
            .help = "Bad pixels mask, HDF5 file"
        }

        output {
          events = None
            .type = path
            .help = "The output HDF5 file"
        }
    """
)


def run():
    # Parse command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--debug", action="store_const", const=True)
    parser.add_argument("phil", nargs="+", action=CheckFileExt)

    args = parser.parse_args()

    cl = phil_scope.command_line_argument_interpreter()
    working_phil = phil_scope.fetch(cl.process_and_fetch(args.phil))
    params = working_phil.extract()

    images_to_events(params)


if __name__ == "__main__":
    with show_mail_on_error():
        run()

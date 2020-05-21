# LIBTBX_SET_DISPATCHER_NAME dev.dials.generate_events
# DIALS_ENABLE_COMMAND_LINE_COMPLETION

import argparse

from libtbx import phil
from dials.util import show_mail_on_error
from dials.algorithms.event_mode import image_to_events

"""

"""

# Create phil parameters
phil_scope = phil.parse(
    """
        input {
          experiments = None
            .type = path
            .multiple = True
            .help = "Input experiment files"

          mask = None
          `.type = path
           .multiple = True
           .help = "Bad pixels mask, HDF5 file"

          image_range = None
            .type = ints(value_min = 0, size = 2)
            .help = "Image range e.g. 1,1000"
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
    parser.add_argument("phil", nargs="+")

    args = parser.parse_args()
    # debug = args.debug

    cl = phil_scope.command_line_argument_interpreter()
    working_phil = phil_scope.fetch(cl.process_and_fetch(args.phil))
    params = working_phil.extract()

    image_to_events.images_to_events(params).run()


if __name__ == "__main__":
    with show_mail_on_error():
        run()

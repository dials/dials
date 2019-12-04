from __future__ import absolute_import, division, print_function

from libtbx.phil import parse
from dials.util import Sorry, show_mail_on_error
from dials.util.options import OptionParser, flatten_reflections, flatten_experiments
from dials.array_family import flex
from dxtbx.model.experiment_list import ExperimentList

from dials.util.multi_dataset_handling import assign_unique_identifiers

help_message = """

Utility to wrangle experiments from indexing in preparation for refinement - 
will look to see that the sample is centred, that there are observations for 
the full range of every scan and rewrite the experiment lists if not. 

  dials.bullseye indexed.expt indexed.refl

"""


phil_scope = parse(
    """
output {
  experiments = wrangled.expt
    .type = path
    .help = "Filename for fixed experiment list"

  reflections = wrangled.refl
    .type = path
    .help = "Filename for fixed reflection list"
}
"""
)


def reflections_fill_scan(reflections, scan, step_degrees):
    """Verify that the reflections defined in the input have Z centroids which 
    fill the scan."""

    z = reflections["xyzobs.px.value"].parts()[2]

    i0, i1 = scan.get_array_range()
    osc = scan.get_oscillation()[1]

    step = int(round(step_degrees / osc))

    for i in range(i0, i1, step):
        if ((z >= i) & (z < (i + step))).count(True) == 0:
            return False

    return True


class Protocol(object):
    def __init__(self):
        usage = "usage: dials.bullseye [options] indexed.expt indexed.refl"

        self.parser = OptionParser(
            usage=usage,
            phil=phil_scope,
            read_reflections=True,
            read_experiments=True,
            check_format=False,
            epilog=help_message,
        )

    def run(self):
        params, _ = self.parser.parse_args(show_diff_phil=True)

        if not params.input.experiments:
            print("No experiments found in the input")
            self.parser.print_help()
            return
        if not params.input.reflections:
            print("No reflections found in the input")
            self.parser.print_help()
            return

        if len(params.input.reflections) != len(params.input.experiments):
            print("Cardinality error")
            self.parser.print_help()
            return

        experiments = flatten_experiments(params.input.experiments)
        reflections = flatten_reflections(params.input.reflections)[0]

        for j, e in enumerate(experiments):
            sel = reflections.select(reflections["id"] == j)
            i0, i1 = e.scan.get_image_range()
            print("Experiment %d has %d reflections" % (j, sel.size()))
            print(
                "Fill scan from %d to %d: %s"
                % (i0, i1, reflections_fill_scan(sel, e.scan, 5.0))
            )


if __name__ == "__main__":
    with show_mail_on_error():
        protocol = Protocol()
        protocol.run()

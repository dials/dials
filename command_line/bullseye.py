from __future__ import absolute_import, division, print_function

from libtbx.phil import parse
from dials.util import Sorry, show_mail_on_error
from dials.util.options import OptionParser, flatten_reflections, flatten_experiments
from dials.array_family import flex
from dxtbx.model.experiment_list import ExperimentList


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


class Script(object):
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
        reflections = flatten_reflections(params.input.reflections)


if __name__ == "__main__":
    with show_mail_on_error():
        script = Script()
        script.run()

from __future__ import absolute_import, division, print_function

from itertools import groupby
import copy
from libtbx.phil import parse
from dials.util import show_mail_on_error
from dials.util.options import OptionParser, flatten_reflections, flatten_experiments
from dials.array_family import flex

help_message = """

Utility to wrangle experiments from indexing in preparation for refinement -
will look to see that the sample is centred, that there are observations for
the full range of every scan and rewrite the experiment lists if not.

  dials.expeditor indexed.expt indexed.refl

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


def select_scans_from_reflections(reflections, scan):
    """Determine a list of valid scans where reflection shoeboxes are seen,
    from within the bounds of the input scan."""

    bbox = reflections["bbox"].parts()

    z0, z1 = bbox[4], bbox[5]
    i0, i1 = scan.get_array_range()

    coverage = flex.int(i1 - i0, 0)

    for _z0, _z1 in zip(z0, z1):
        for j in range(_z0, _z1):
            coverage[j] += 1

    filled = (coverage > 0).iselection()

    scans = []

    for k, g in groupby(enumerate(filled), lambda n: n[0] - n[1]):
        l = list(g)
        s = copy.deepcopy(scan)
        s.set_image_range((l[0][1] + 1, l[-1][1] + 1))
        scans.append(s)

    return scans


def reflections_fill_scan(reflections, scan, step_degrees):
    """Verify that the reflections defined in the input have Z centroids which
    fill the scan."""

    z = reflections["xyzobs.px.value"].parts()[2]

    i0, i1 = scan.get_array_range()
    slots = int(round((i1 - i0) * scan.get_oscillation()[1] / step_degrees))

    fill = flex.histogram(z, data_min=i0, data_max=i1, n_slots=slots)

    for s in fill.slots():
        if s == 0:
            return False

    return True


class Protocol(object):
    def __init__(self):
        usage = "usage: dials.expeditor [options] indexed.expt indexed.refl"

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
            scans = select_scans_from_reflections(sel, e.scan)
            for s in scans:
                print("Image range: %d %d" % s.get_image_range())


if __name__ == "__main__":
    with show_mail_on_error():
        protocol = Protocol()
        protocol.run()

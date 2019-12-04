from __future__ import absolute_import, division, print_function

from itertools import groupby
from tabulate import tabulate
import copy
from libtbx.phil import parse
from dials.util import show_mail_on_error
from dials.util.options import OptionParser, flatten_reflections, flatten_experiments
from dials.array_family import flex
from dxtbx.model.experiment_list import ExperimentList

help_message = """

Utility to wrangle experiments from indexing in preparation for refinement -
will look to see that the sample is centred, that there are observations for
the full range of every scan and rewrite the experiment lists if not.

  dials.expeditor indexed.expt indexed.refl

"""


phil_scope = parse(
    """
output {
  experiments = expedited.expt
    .type = path
    .help = "Filename for fixed experiment list"

  reflections = expedited.refl
    .type = path
    .help = "Filename for fixed reflection list"
}
"""
)


def model_connectivity(experiments):
    def model_connectivity_impl(experiments, model):
        text = [""]
        text.append("%s:" % model.capitalize())
        models = getattr(experiments, "%ss" % model)()
        rows = [[""] + [str(j) for j in range(len(models))]]
        for j, e in enumerate(experiments):
            row = ["Experiment %d" % j]
            for m in models:
                if getattr(e, model) is m:
                    row.append("x")
                else:
                    row.append(".")
            rows.append(row)
        text.append(tabulate(rows, tablefmt="plain"))
        return text

    if len(experiments) == 1:
        return ""

    text = []
    text.append("Experiment / Models")
    text.extend(model_connectivity_impl(experiments, "detector"))
    text.extend(model_connectivity_impl(experiments, "crystal"))
    text.extend(model_connectivity_impl(experiments, "beam"))

    return "\n".join(text)


def select_scans_from_reflections(reflections, scan):
    """Determine a list of valid scans where reflection shoeboxes are seen,
    from within the bounds of the input scan."""

    bbox = reflections["bbox"].parts()

    z0, z1 = bbox[4], bbox[5]
    i0, i1 = scan.get_array_range()
    o0, o1 = scan.get_oscillation()

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
        s.set_oscillation((o0 + o1 * (l[0][1] + 1 - i0), o1))
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

        print(model_connectivity(experiments))
        print()

        z = reflections["xyzobs.px.value"].parts()[2]

        experiments_out = ExperimentList()
        reflections_id = flex.int(reflections.size(), -424242)
        reflections_id.set_selected(reflections["id"] == -1, -1)

        crystal_scan = {}

        for j, e in enumerate(experiments):
            if e.crystal in crystal_scan:
                if e.scan is not crystal_scan[e.crystal]:
                    e.crystal = copy.deepcopy(e.crystal)
            else:
                crystal_scan[e.crystal] = e.scan

            sel = reflections.select(reflections["id"] == j)
            i0, i1 = e.scan.get_image_range()
            nref = sel.size()
            scans = select_scans_from_reflections(sel, e.scan)

            assert len(scans) > 0

            eid = len(experiments_out)

            e.scan = scans[0]
            experiments_out.append(e)

            # rewrite id
            i0, i1 = e.scan.get_array_range()
            sel = (reflections["id"] == j) & (z >= i0) & (z <= i1)
            reflections_id.set_selected(sel, eid)

            print(
                "Output experiment %d has %d / %d reflections"
                % (eid, sel.count(True), nref)
            )

            for k, s in enumerate(scans[1:]):
                f = copy.deepcopy(e)
                f.scan = s
                i0, i1 = s.get_array_range()
                sel = (reflections["id"] == j) & (z >= i0) & (z <= i1)
                print(
                    "Output experiment %d has %d / %d reflections"
                    % (eid + 1 + k, sel.count(True), nref)
                )
                reflections_id.set_selected(sel, eid + 1 + k)
                experiments_out.append(f)

        assert reflections_id.count(-424242) == 0

        reflections["id"] = reflections_id

        print()
        print(model_connectivity(experiments_out))
        experiments_out.as_file(params.output.experiments)
        reflections.as_file(params.output.reflections)


if __name__ == "__main__":
    with show_mail_on_error():
        protocol = Protocol()
        protocol.run()

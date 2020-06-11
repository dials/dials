# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# DIALS_ENABLE_COMMAND_LINE_COMPLETION
from __future__ import absolute_import, division, print_function

import copy
import sys

from scitbx.array_family import flex
import wxtbx.app

import dials.util.log
from dials.util.reciprocal_lattice.viewer import ReciprocalLatticeViewer, phil_scope
from dials.util.options import OptionParser, reflections_and_experiments_from_files

help_message = """
Visualise the strong spots from spotfinding in reciprocal space.

Examples::

  dials.reciprocal_lattice_viewer imported.expt strong.refl

  dials.reciprocal_lattice_viewer indexed.expt indexed.refl
"""


def run(args):
    usage = "dials.reciprocal_lattice_viewer [options] models.expt observations.refl"

    parser = OptionParser(
        usage=usage,
        phil=phil_scope,
        read_experiments=True,
        read_reflections=True,
        check_format=False,
        epilog=help_message,
    )

    params, options = parser.parse_args(show_diff_phil=True)
    reflections, experiments = reflections_and_experiments_from_files(
        params.input.reflections, params.input.experiments
    )

    if len(experiments) == 0 or len(reflections) == 0:
        parser.print_help()
        exit(0)

    if len(reflections) > 1:
        assert len(reflections) == len(experiments)
        for i in range(len(reflections)):
            reflections[i]["imageset_id"] = flex.int(len(reflections[i]), i)
            if i > 0:
                reflections[0].extend(reflections[i])
    elif "imageset_id" not in reflections[0]:
        reflections[0]["imageset_id"] = reflections[0]["id"]

    reflections = reflections[0]

    a = wxtbx.app.CCTBXApp(0)
    a.settings = params
    f = ReciprocalLatticeViewer(
        None,
        -1,
        "Reflection data viewer",
        size=(1024, 768),
        settings=copy.deepcopy(params),
    )
    f.load_models(experiments, reflections)
    f.Show()
    a.SetTopWindow(f)
    a.MainLoop()


if __name__ == "__main__":
    dials.util.log.print_banner()
    run(sys.argv[1:])

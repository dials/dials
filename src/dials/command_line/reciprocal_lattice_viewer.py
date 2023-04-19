# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# DIALS_ENABLE_COMMAND_LINE_COMPLETION
# LIBTBX_SET_DISPATCHER_NAME dials.reciprocal_lattice_viewer
# LIBTBX_SET_DISPATCHER_NAME dials.rlv

from __future__ import annotations

import copy
import os

import wxtbx.app
from scitbx.array_family import flex

import dials.util.log
from dials.util.options import ArgumentParser, reflections_and_experiments_from_files
from dials.util.reciprocal_lattice.viewer import ReciprocalLatticeViewer, phil_scope

help_message = """
Visualise the strong spots from spotfinding in reciprocal space.

Examples::

  dials.reciprocal_lattice_viewer imported.expt strong.refl

  dials.reciprocal_lattice_viewer indexed.expt indexed.refl
"""


@dials.util.show_mail_handle_errors()
def run(args=None):
    dials.util.log.print_banner()
    usage = "dials.reciprocal_lattice_viewer [options] models.expt observations.refl"

    parser = ArgumentParser(
        usage=usage,
        phil=phil_scope,
        read_experiments=True,
        read_reflections=True,
        check_format=False,
        epilog=help_message,
    )

    params, options = parser.parse_args(args, show_diff_phil=True)
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
        os.path.realpath(params.input.reflections[0].filename),
        size=(1024, 768),
        settings=copy.deepcopy(params),
    )
    f.load_models(experiments, reflections)
    f.Show()
    a.SetTopWindow(f)
    a.MainLoop()


if __name__ == "__main__":
    run()

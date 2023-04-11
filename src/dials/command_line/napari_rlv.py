# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_SET_DISPATCHER_NAME dev.dials.napari_rlv

from __future__ import annotations

import copy
import os
import sys

from packaging import version

from scitbx.array_family import flex

import dials.util.log
from dials.util.napari_rlv.viewer import ReciprocalLatticeViewer, phil_scope
from dials.util.options import ArgumentParser, reflections_and_experiments_from_files

try:
    import napari
except ImportError:
    napari = None

help_message = """
Visualise the strong spots from spotfinding in reciprocal space.

Examples::

  dev.dials.napari_rlv imported.expt strong.refl

  dev.dials.napari_rlv indexed.expt indexed.refl
"""


@dials.util.show_mail_handle_errors()
def run(args=None):
    dials.util.log.print_banner()
    usage = "dev.dials.napari_rlv [options] models.expt observations.refl"

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

    napari_viewer = napari.Viewer(
        title=os.path.realpath(params.input.reflections[0].filename),
        ndisplay=3,
    )
    rlv = ReciprocalLatticeViewer(
        None,
        -1,
        os.path.realpath(params.input.reflections[0].filename),
        size=(1024, 768),
        settings=copy.deepcopy(params),
        napari_viewer=napari_viewer,
    )
    rlv.load_models(experiments, reflections)
    rlv.add_rlv_widgets()
    rlv.update_layers()

    # Set rotation around the origin
    napari_viewer.camera.center = (0, 0, 0)

    # Make the last-added relps layer active rather than the rotation axis
    napari_viewer.layers.selection.active = napari_viewer.layers[0]

    napari.run()


if __name__ == "__main__":
    if not napari:
        sys.exit("Please install napari")
    if version.parse(napari.__version__) <= version.parse("0.4.15"):
        sys.exit("Please install napari >= 0.4.16")

    run()

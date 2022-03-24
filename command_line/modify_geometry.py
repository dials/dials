from __future__ import annotations

from typing import List

import libtbx.phil
from dxtbx.model import ExperimentList

import dials.util
from dials.command_line.dials_import import ManualGeometryUpdater
from dials.util.options import ArgumentParser, flatten_experiments

help_message = """
"""

phil_scope = libtbx.phil.parse(
    """
include scope dials.util.options.geometry_phil_scope
output {
  experiments = modified.expt
    .type = path
}
""",
    process_includes=True,
)


def update(
    experiments: ExperimentList, new_params: libtbx.phil.scope_extract
) -> ExperimentList:

    """
    Modify detector, beam, goniometer and scan in experiments with the values in new_params
    """

    update_geometry = ManualGeometryUpdater(new_params)

    imagesets = experiments.imagesets()

    for imageset in imagesets:
        imageset_new = update_geometry(imageset)
        imageset.set_detector(imageset_new.get_detector())
        imageset.set_beam(imageset_new.get_beam())
        imageset.set_goniometer(imageset_new.get_goniometer())
        imageset.set_scan(imageset_new.get_scan())

    return experiments


@dials.util.show_mail_handle_errors()
def run(args: List[str] = None, phil: libtbx.phil.scope = phil_scope) -> None:

    usage = "dials.modify_geometry [options] models.expt"

    parser = ArgumentParser(
        usage=usage,
        phil=phil_scope,
        read_experiments=True,
        check_format=False,
        epilog=help_message,
    )

    params, _ = parser.parse_args(args, show_diff_phil=True)

    experiments = flatten_experiments(params.input.experiments)

    if len(experiments) == 0:
        parser.print_help()
        exit(0)

    new_experiments = update(experiments, params)

    if len(new_experiments):
        print(f"Saving modified experiments to {params.output.experiments}")
        new_experiments.as_file(params.output.experiments)


if __name__ == "__main__":
    run()

from __future__ import annotations

import libtbx.phil
from dxtbx.model import ExperimentList
from dxtbx.model.crystal import CrystalFactory

import dials.util
from dials.command_line.dials_import import ManualGeometryUpdater
from dials.util.options import ArgumentParser, flatten_experiments

help_message = """
"""

phil_scope = libtbx.phil.parse(
    """
include scope dials.util.options.geometry_phil_scope
include scope dxtbx.model.crystal.crystal_phil_scope
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
    Modify the models in experiments with the values in new_params
    """

    update_geometry = ManualGeometryUpdater(new_params)

    for experiment in experiments:
        imageset = update_geometry(experiment.imageset)
        experiment.imageset = imageset
        experiment.detector = imageset.get_detector()
        experiment.beam = imageset.get_beam()
        experiment.goniometer = imageset.get_goniometer()
        experiment.scan = imageset.get_scan()
        experiment.scan.set_valid_image_ranges(experiment.identifier, [])
        crystal = CrystalFactory.from_phil(new_params, experiment.crystal)
        experiment.crystal = crystal
    return experiments


@dials.util.show_mail_handle_errors()
def run(args: list[str] = None, phil: libtbx.phil.scope = phil_scope) -> None:
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

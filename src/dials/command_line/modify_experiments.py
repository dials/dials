from __future__ import annotations

from copy import deepcopy

import libtbx.phil
from dxtbx.model import ExperimentList
from dxtbx.model.crystal import CrystalFactory

import dials.util
from dials.command_line.dials_import import ManualGeometryUpdater
from dials.util.options import ArgumentParser, flatten_experiments

help_message = """

This program modifies the experimental models of the experiments in an experiment list.

Examples::

  dials.modify_experiments models.expt distance=100.0

  dials.modify_experiments models.expt select_experiments=0,1 A_matrix=-0.076948,0.058256,0.104294,-0.010462,0.113451,-0.081650,-0.112936,-0.050201,-0.063496
"""

phil_scope = libtbx.phil.parse(
    """
include scope dials.util.options.geometry_phil_scope
include scope dxtbx.model.crystal.crystal_phil_scope
select_experiments = None
    .type = ints
    .help = "A list of experiment ids to select for modification. If None, all experiments are modified."
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
    if new_params.select_experiments is None:
        new_params.select_experiments = list(range(len(experiments)))
    else:
        # Copy the selected experiments to ensure they do not share models with
        # the original experiments
        for iexp in new_params.select_experiments:
            if iexp < 0 or iexp >= len(experiments):
                raise ValueError(
                    f"Experiment index {iexp} is out of range for the number of experiments ({len(experiments)})"
                )
            experiments[iexp] = deepcopy(experiments[iexp])
    for iexp in new_params.select_experiments:
        experiment = experiments[iexp]
        imageset = update_geometry(experiment.imageset)
        experiment.imageset = imageset
        experiment.detector = imageset.get_detector()
        experiment.beam = imageset.get_beam()
        experiment.goniometer = imageset.get_goniometer()
        experiment.scan = imageset.get_scan()
        experiment.scan.set_valid_image_ranges(experiment.identifier, [])
        crystal = CrystalFactory.from_phil(new_params, experiment.crystal)
        experiment.crystal = crystal
        experiments[iexp] = experiment
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

from __future__ import annotations

import copy
import warnings

import libtbx.phil
from dxtbx.model import ExperimentList
from dxtbx.serialize import load

import dials.util
from dials.command_line.dials_import import ManualGeometryUpdater
from dials.util.options import ArgumentParser, flatten_experiments

help_message = """
"""

phil_scope = libtbx.phil.parse(
    """
include scope dials.util.options.geometry_phil_scope
input {
   reference_geometry = None
       .type = path
        help = "Experimental geometry from this models.expt "
              "will override the geometry from the input file"
    use_beam_reference = True
      .type = bool
      .expert_level = 2
      .help = "If True, the beam from reference_geometry will override "
              "the beam from the input file."
    use_gonio_reference = True
      .type = bool
      .expert_level = 2
      .help = "If True, the goniometer from reference_geometry will override "
              "the goniometer from the input file."
    use_detector_reference = True
      .type = bool
      .expert_level = 2
      .help = "If True, the detector from reference_geometry will override "
              "the detector from the input file"
}
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

    if params.input.reference_geometry:
        ref_expts = load.experiment_list(params.input.reference_geometry)
        if params.input.use_detector_reference:
            ref_det = ref_expts[0].detector
            current_detectors = experiments.detectors()
            # we want to retain the structure of the models i.e. shared or not.
            new_detectors = [
                copy.deepcopy(ref_det) for _ in range(len(current_detectors))
            ]
            for expt in experiments:
                expt.detector = new_detectors[current_detectors.index(expt.detector)]
        if params.input.use_beam_reference:
            ref_beam = ref_expts[0].beam
            current_beams = experiments.beams()
            # we want to retain the structure of the models i.e. shared or not.
            new_beams = [copy.deepcopy(ref_beam) for _ in range(len(current_beams))]
            for expt in experiments:
                expt.beam = new_beams[current_beams.index(expt.beam)]
        if params.input.use_gonio_reference:
            ref_gonio = ref_expts[0].goniometer
            current_gonios = experiments.goniometers()
            # we want to retain the structure of the models i.e. shared or not.
            new_gonios = [copy.deepcopy(ref_gonio) for _ in range(len(current_gonios))]
            for expt in experiments:
                expt.goniometer = new_gonios[current_gonios.index(expt.goniometer)]

    # update with any manual parameters set.
    new_experiments = update(experiments, params)

    if len(new_experiments):
        print(f"Saving modified experiments to {params.output.experiments}")
        new_experiments.as_file(params.output.experiments)


if __name__ == "__main__":
    warnings.warn(
        "dials.modify_geometry is deprecated, please use dials.modify_experiments instead.\n",
        DeprecationWarning,
        stacklevel=1,
    )
    run()

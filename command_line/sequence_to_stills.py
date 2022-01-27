"""
Split a sequence into a set of stills.

Example:

  dials.sequence_to_stills sequence.expt sequence.refl
"""


from __future__ import annotations

import logging

from dxtbx.model import MosaicCrystalSauter2014
from dxtbx.model.experiment_list import Experiment, ExperimentList
from libtbx.phil import parse
from scitbx import matrix

from dials.algorithms.refinement.prediction.managed_predictors import (
    ExperimentsPredictorFactory,
)
from dials.array_family import flex
from dials.model.data import Shoebox
from dials.util import show_mail_handle_errors
from dials.util.options import ArgumentParser, reflections_and_experiments_from_files

logger = logging.getLogger("dials.command_line.sequence_to_stills")

# The phil scope
phil_scope = parse(
    """
output {
  experiments = stills.expt
    .type = str
    .help = "Filename for the experimental models that have been converted to stills"
  reflections = stills.refl
    .type = str
    .help = "Filename for the reflection tables with split shoeboxes (3D to 2D)"
  domain_size_ang = None
    .type = float
    .help = "Override for domain size. If None, use the crystal's domain size, if"
            "available"
  half_mosaicity_deg = None
    .type = float
    .help = "Override for mosaic angle. If None, use the crystal's mosaic angle, if"
            "available"
}
max_scan_points = None
  .type = int
  .expert_level = 2
  .help = Limit number of scan points
"""
)


def sequence_to_stills(experiments, reflections, params):
    assert len(reflections) == 1
    reflections = reflections[0]

    new_experiments = ExperimentList()
    new_reflections = flex.reflection_table()

    # This is the subset needed to integrate
    for key in [
        "id",
        "imageset_id",
        "shoebox",
        "bbox",
        "intensity.sum.value",
        "intensity.sum.variance",
        "entering",
        "flags",
        "miller_index",
        "panel",
        "xyzobs.px.value",
        "xyzobs.px.variance",
    ]:
        if key in reflections:
            new_reflections[key] = type(reflections[key])()
        elif key == "imageset_id":
            assert len(experiments.imagesets()) == 1
            reflections["imageset_id"] = flex.int(len(reflections), 0)
            new_reflections["imageset_id"] = flex.int()
        elif key == "entering":
            reflections["entering"] = flex.bool(len(reflections), False)
            new_reflections["entering"] = flex.bool()
        else:
            raise RuntimeError(f"Expected key not found in reflection table: {key}")

    for expt_id, experiment in enumerate(experiments):
        # Get the goniometer setting matrix
        goniometer_setting_matrix = matrix.sqr(
            experiment.goniometer.get_setting_rotation()
        )
        goniometer_axis = matrix.col(experiment.goniometer.get_rotation_axis())
        step = experiment.scan.get_oscillation()[1]

        refls = reflections.select(reflections["id"] == expt_id)
        _, _, _, _, z1, z2 = refls["bbox"].parts()

        # Create an experiment for each scanpoint.
        # Note that a simplification of the use of scan-points here introduces a
        # (small) error. The scan-varying crystal model has 1 more scan-point
        # than the scan has images because scan-points are taken at the boundaries
        # between images, including the scan extrema. This code assumes that
        # the crystal model at the start of each image applies to the whole
        # image and ignores the final scan-point.
        start, stop = experiment.scan.get_array_range()
        for i_array in range(start, stop):
            if params.max_scan_points and i_array >= params.max_scan_points:
                break
            # Shift array position to scan-point index
            i_scan_point = i_array - start

            # The A matrix is the goniometer setting matrix for this scan point
            # times the scan varying A matrix at this scan point. Note, the
            # goniometer setting matrix for scan point zero will be the identity
            # matrix and represents the beginning of the oscillation.
            # For stills, the A matrix needs to be positioned in the midpoint of an
            # oscillation step. Hence, here the goniometer setting matrix is rotated
            # by a further half oscillation step.
            A = (
                goniometer_axis.axis_and_angle_as_r3_rotation_matrix(
                    angle=experiment.scan.get_angle_from_array_index(i_array)
                    + (step / 2),
                    deg=True,
                )
                * goniometer_setting_matrix
                * matrix.sqr(experiment.crystal.get_A_at_scan_point(i_array))
            )
            crystal = MosaicCrystalSauter2014(experiment.crystal)
            crystal.set_A(A)

            # Copy in mosaic parameters if available
            if params.output.domain_size_ang is None and hasattr(
                experiment.crystal, "get_domain_size_ang"
            ):
                crystal.set_domain_size_ang(experiment.crystal.get_domain_size_ang())
            elif params.output.domain_size_ang is not None:
                crystal.set_domain_size_ang(params.output.domain_size_ang)

            if params.output.half_mosaicity_deg is None and hasattr(
                experiment.crystal, "get_half_mosaicity_deg"
            ):
                crystal.set_half_mosaicity_deg(
                    experiment.crystal.get_half_mosaicity_deg()
                )
            elif params.output.half_mosaicity_deg is not None:
                crystal.set_half_mosaicity_deg(params.output.half_mosaicity_deg)

            new_experiment = Experiment(
                detector=experiment.detector,
                beam=experiment.beam,
                crystal=crystal,
                imageset=experiment.imageset.as_imageset()[
                    i_scan_point : i_scan_point + 1
                ],
            )
            new_experiments.append(new_experiment)

            # Each reflection in a 3D shoebox can be found on multiple images.
            # Slice the reflections such that any reflection on this scan point
            # is included with this image
            new_id = len(new_experiments) - 1
            subrefls = refls.select((i_scan_point >= z1) & (i_scan_point < z2))
            for refl in subrefls.rows():
                assert i_scan_point in range(*refl["bbox"][4:6])

                new_sb = Shoebox()
                start = i_scan_point - refl["bbox"][4]  # z1
                new_sb.data = refl["shoebox"].data[start : start + 1, :, :]
                new_sb.background = refl["shoebox"].background[start : start + 1, :, :]
                new_sb.mask = refl["shoebox"].mask[start : start + 1, :, :]
                intensity = new_sb.summed_intensity()
                new_sb.bbox = tuple(
                    list(refl["bbox"])[0:4] + [0, 1]
                )  # keep the original shoebox but reset the z values
                new_sb.panel = refl["panel"]
                new_refl = {}
                new_refl["id"] = new_refl["imageset_id"] = new_id
                new_refl["shoebox"] = new_sb
                new_refl["bbox"] = new_sb.bbox
                new_refl["intensity.sum.value"] = intensity.observed.value
                new_refl["intensity.sum.variance"] = intensity.observed.variance
                for key in ["entering", "flags", "miller_index", "panel"]:
                    new_refl[key] = refl[key]
                centroid = new_sb.centroid_foreground_minus_background()
                new_refl["xyzobs.px.value"] = centroid.px.position
                new_refl["xyzobs.px.variance"] = centroid.px.variance
                new_reflections.append({})
                for key in new_refl:
                    new_reflections[key][-1] = new_refl[key]

    # Re-predict using the reflection slices and the stills predictors
    ref_predictor = ExperimentsPredictorFactory.from_experiments(
        new_experiments, force_stills=new_experiments.all_stills()
    )
    new_reflections = ref_predictor(new_reflections)

    return (new_experiments, new_reflections)


@show_mail_handle_errors()
def run(args=None, phil=phil_scope):
    """
    Validate the arguments and load experiments/reflections for sequence_to_stills

    Arguments:
        args: The command line arguments to use. Defaults to sys.argv[1:]
        phil: The phil_scope. Defaults to the master phil_scope for this program
    """
    # The script usage
    usage = "usage: dials.sequence_to_stills [options] [param.phil] models.expt reflections.refl"

    # Create the parser
    parser = ArgumentParser(
        usage=usage,
        phil=phil,
        read_experiments=True,
        read_reflections=True,
        check_format=False,
        epilog=__doc__,
    )
    params, options = parser.parse_args(args=args, show_diff_phil=True)

    # Try to load the models and data
    if not params.input.experiments or not params.input.reflections:
        parser.print_help()
        return

    reflections, experiments = reflections_and_experiments_from_files(
        params.input.reflections, params.input.experiments
    )

    (new_experiments, new_reflections) = sequence_to_stills(
        experiments, reflections, params
    )
    # Write out the output experiments, reflections
    new_experiments.as_file(params.output.experiments)
    new_reflections.as_file(params.output.reflections)


if __name__ == "__main__":
    run()

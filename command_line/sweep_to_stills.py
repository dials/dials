#!/usr/bin/env python
#
# LIBTBX_SET_DISPATCHER_NAME dials.sweep_to_stills

from __future__ import absolute_import, division, print_function

import logging

from dxtbx.model.experiment_list import Experiment, ExperimentList
from scitbx import matrix
from dxtbx.model import MosaicCrystalSauter2014
from dials.model.data import Shoebox
from dials.algorithms.refinement.prediction.managed_predictors import (
    ExperimentsPredictorFactory,
)
from dials.util import show_mail_on_error
from libtbx.phil import parse

logger = logging.getLogger("dials.command_line.sweep_to_stills")

help_message = """

Split a sweep into a set of stills.

Example::

  dials.sweep_to_stills sweep.expt sweep.refl

"""

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
"""
)


class Script(object):
    """A class for running the script."""

    def __init__(self):
        """Initialise the script."""
        from dials.util.options import OptionParser
        import libtbx.load_env

        # The script usage
        usage = (
            "usage: %s [options] [param.phil] filenames" % libtbx.env.dispatcher_name
        )

        # Create the parser
        self.parser = OptionParser(
            usage=usage,
            phil=phil_scope,
            read_experiments=True,
            read_reflections=True,
            check_format=False,
            epilog=help_message,
        )

    def run(self):
        """Execute the script."""
        from dials.util.options import flatten_experiments, flatten_reflections
        from dials.array_family import flex

        # Parse the command line
        params, options = self.parser.parse_args(show_diff_phil=True)

        # Try to load the models and data
        if not params.input.experiments:
            print("No Experiments found in the input")
            self.parser.print_help()
            return
        if not params.input.reflections:
            print("No reflections found in the input")
            self.parser.print_help()
            return

        experiments = flatten_experiments(params.input.experiments)
        reflections = flatten_reflections(params.input.reflections)
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
            new_reflections[key] = type(reflections[key])()

        for expt_id, experiment in enumerate(experiments):
            # Get the goniometr setting matrix
            goniometer_setting_matrix = matrix.sqr(
                experiment.goniometer.get_setting_rotation()
            )
            goniometer_axis = matrix.col(experiment.goniometer.get_rotation_axis())
            step = experiment.scan.get_oscillation()[1]

            refls = reflections.select(reflections["id"] == expt_id)
            _, _, _, _, z1, z2 = refls["bbox"].parts()

            # Create an experiment for each scanpoint
            for i_scan_point in range(*experiment.scan.get_array_range()):
                # The A matrix is the goniometer setting matrix for this scan point
                # times the scan varying A matrix at this scan point. Note, the
                # goniometer setting matrix for scan point zero will be the identity
                # matrix and represents the beginning of the oscillation.
                # For stills, the A matrix needs to be positioned in the midpoint of an
                # oscillation step. Hence, here the goniometer setting matrixis rotated
                # by a further half oscillation step.
                A = (
                    goniometer_axis.axis_and_angle_as_r3_rotation_matrix(
                        angle=experiment.scan.get_angle_from_array_index(i_scan_point)
                        + (step / 2),
                        deg=True,
                    )
                    * goniometer_setting_matrix
                    * matrix.sqr(experiment.crystal.get_A_at_scan_point(i_scan_point))
                )
                crystal = MosaicCrystalSauter2014(experiment.crystal)
                crystal.set_A(A)

                # Copy in mosaic parameters if available
                if params.output.domain_size_ang is None and hasattr(
                    experiment.crystal, "get_domain_size_ang"
                ):
                    crystal.set_domain_size_ang(
                        experiment.crystal.get_domain_size_ang()
                    )
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
                    new_sb.background = refl["shoebox"].background[
                        start : start + 1, :, :
                    ]
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

            new_experiments.as_file(params.output.experiments)
            new_reflections.as_file(params.output.reflections)


if __name__ == "__main__":
    with show_mail_on_error():
        script = Script()
        script.run()

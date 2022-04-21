# LIBTBX_SET_DISPATCHER_NAME dials.complete_full_sphere


from __future__ import annotations

import logging
import math
import sys

from cctbx import crystal, miller
from dxtbx.model import ExperimentList
from libtbx.phil import parse
from scitbx import matrix

from dials.algorithms.refinement import rotation_decomposition
from dials.algorithms.shadowing.filter import filter_shadowed_reflections
from dials.array_family import flex
from dials.util import log, show_mail_handle_errors
from dials.util.options import ArgumentParser, flatten_experiments

logger = logging.getLogger("dials.command_line.complete_full_sphere")

help_message = """
This program attempts to compute a sample realignment to measure the blind
region of reciprocal space, given an already recorded data set.

dials.complete_full_sphere [resolution=1.6] models.expt
"""

phil_scope = parse(
    """
resolution = 0.0
  .type = float
  .help = "Resolution of diffraction for blind region calculation"
shadow = True
  .type = bool
  .help = "Consider shadowing in calculating overall completeness"
"""
)


class Script:
    """A class for running the script."""

    def __init__(self):
        """Initialise the script."""
        # The script usage
        usage = "usage: dials.complete_full_sphere [options] "

        # Create the parser
        self.parser = ArgumentParser(
            usage=usage,
            phil=phil_scope,
            epilog=help_message,
            check_format=True,
            read_experiments=True,
        )

    def run(self, args=None):
        params, options = self.parser.parse_args(args, show_diff_phil=True)
        log.config(logfile="dials.complete_full_sphere.log")

        model_shadow = params.shadow

        experiments = flatten_experiments(params.input.experiments)

        if len(experiments) != 1:
            self.parser.print_help()
            return

        expt = experiments[0]

        axes = expt.goniometer.get_axes()

        if len(axes) != 3:
            sys.exit("This will only work with 3-axis goniometers")

        if not expt.imageset.reader().get_format():
            sys.exit("This will only work with images available")

        if not expt.imageset.reader().get_format().get_goniometer_shadow_masker():
            model_shadow = False

        beam = expt.beam
        det = expt.detector

        if params.resolution:
            resolution = params.resolution
        else:
            resolution = det.get_max_inscribed_resolution(expt.beam.get_s0())

        # at this point, predict all of the reflections in the scan possible (i.e.
        # extend scan to 360 degrees) - this points back to expt

        self.make_scan_360(expt.scan)

        # now get a full set of all unique miller indices
        all_indices = miller.build_set(
            crystal_symmetry=crystal.symmetry(
                space_group=expt.crystal.get_space_group(),
                unit_cell=expt.crystal.get_unit_cell(),
            ),
            anomalous_flag=True,
            d_min=resolution,
        )

        if model_shadow:
            obs, shadow = self.predict_to_miller_set_with_shadow(expt, resolution)
        else:
            obs = self.predict_to_miller_set(expt, resolution)

        logger.info(
            "Fraction of unique observations at datum: %.1f%%",
            100.0 * len(obs.indices()) / len(all_indices.indices()),
        )

        missing = all_indices.lone_set(other=obs)

        logger.info("%d unique reflections in blind region", len(missing.indices()))

        e1 = matrix.col(axes[0])
        e2 = matrix.col(axes[1])
        e3 = matrix.col(axes[2])

        s0n = matrix.col(beam.get_s0()).normalize()

        # rotate blind region about beam by +/- two theta
        two_theta = 2.0 * math.asin(0.5 * beam.get_wavelength() / resolution)
        R_ptt = s0n.axis_and_angle_as_r3_rotation_matrix(two_theta)
        R_ntt = s0n.axis_and_angle_as_r3_rotation_matrix(-two_theta)

        # now decompose to e3, e2, e1
        sol_plus = rotation_decomposition.solve_r3_rotation_for_angles_given_axes(
            R_ptt, e3, e2, e1, return_both_solutions=True, deg=True
        )

        sol_minus = rotation_decomposition.solve_r3_rotation_for_angles_given_axes(
            R_ntt, e3, e2, e1, return_both_solutions=True, deg=True
        )

        solutions = []
        if sol_plus:
            solutions.extend(sol_plus)
        if sol_minus:
            solutions.extend(sol_minus)

        if not solutions:
            sys.exit(f"Impossible two theta: {two_theta * 180.0 / math.pi:.3f},")

        logger.info("Maximum two theta: %.3f,", two_theta * 180.0 / math.pi)
        logger.info("%d solutions found", len(solutions))

        names = tuple(
            [n.replace("GON_", "").lower() for n in expt.goniometer.get_names()]
        )

        logger.info(" %8s %8s %8s  coverage expt.expt" % names)
        self.write_expt(experiments, "solution_0.expt")
        for j, s in enumerate(solutions):
            expt.goniometer.set_angles(s)
            if model_shadow:
                obs, shadow = self.predict_to_miller_set_with_shadow(expt, resolution)
            else:
                obs = self.predict_to_miller_set(expt, resolution)
            new = missing.common_set(obs)
            fout = "solution_%d.expt" % (j + 1)
            f = len(new.indices()) / len(missing.indices())

            logger.info("%8.3f %8.3f %8.3f %4.2f %s", s[0], s[1], s[2], f, fout)
            self.write_expt(experiments, fout)

    def make_scan_360(self, scan):
        epochs = scan.get_epochs()
        exposure_times = scan.get_exposure_times()
        image_range = scan.get_image_range()
        oscillation = scan.get_oscillation()

        current = 1 + image_range[1] - image_range[0]
        turn = int(round(360.0 / oscillation[1]))
        extra = turn - current

        for j in range(extra):
            epochs.append(0.0)
            exposure_times.append(0.0)
        image_range = image_range[0], image_range[1] + extra

        scan.set_image_range(image_range)
        scan.set_epochs(epochs)
        scan.set_exposure_times(exposure_times)
        return scan

    def predict_to_miller_set(self, expt, resolution):
        predicted = flex.reflection_table.from_predictions(expt, dmin=resolution)
        hkl = predicted["miller_index"]

        # now get a full set of all unique miller indices
        obs = miller.set(
            crystal_symmetry=crystal.symmetry(
                space_group=expt.crystal.get_space_group(),
                unit_cell=expt.crystal.get_unit_cell(),
            ),
            anomalous_flag=True,
            indices=hkl,
        ).unique_under_symmetry()

        return obs

    def predict_to_miller_set_with_shadow(self, expt, resolution):
        predicted = flex.reflection_table.from_predictions(expt, dmin=resolution)

        # transmogrify this to an ExperimentList from an Experiment
        experiments = ExperimentList()
        experiments.append(expt)
        predicted["id"] = flex.int(predicted.size(), 0)
        shadowed = filter_shadowed_reflections(
            experiments, predicted, experiment_goniometer=True
        )
        predicted = predicted.select(~shadowed)

        hkl = predicted["miller_index"]

        # now get a full set of all unique miller indices
        obs = miller.set(
            crystal_symmetry=crystal.symmetry(
                space_group=expt.crystal.get_space_group(),
                unit_cell=expt.crystal.get_unit_cell(),
            ),
            anomalous_flag=True,
            indices=hkl,
        ).unique_under_symmetry()

        return obs, shadowed

    def write_expt(self, experiments, filename):
        experiments.as_file(filename)


@show_mail_handle_errors()
def run(args=None):
    script = Script()
    script.run(args)


if __name__ == "__main__":
    run()

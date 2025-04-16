"""Managed reflection prediction for refinement.

* ScansRayPredictor adapts DIALS prediction for use in refinement, by keeping
  up to date with the current model geometry
* StillsRayPredictor predicts reflections without a goniometer, under
  the naive assumption that the relp is already in reflecting position
"""

from __future__ import annotations

from math import pi

from dxtbx.model import tof_helpers
from scitbx.array_family import flex

from dials.algorithms.spot_prediction import (
    LaueReflectionPredictor,
    ScanStaticRayPredictor,
)
from dials.algorithms.spot_prediction import ScanStaticReflectionPredictor as sc
from dials.algorithms.spot_prediction import ScanVaryingReflectionPredictor as sv
from dials.algorithms.spot_prediction import StillsReflectionPredictor as st

# constants
TWO_PI = 2.0 * pi


class ScansRayPredictor:
    """
    Predict for a relp based on the current states of models of the
    experimental geometry. This is a wrapper for DIALS' C++
    RayPredictor class, which does the real work. This class keeps track
    of the experimental geometry, and instantiates a RayPredictor when
    required.
    """

    def __init__(self, experiments, sequence_range=(0, 2.0 * pi)):
        """Construct by linking to instances of experimental model classes"""

        self._experiments = experiments
        self._sequence_range = sequence_range

    def __call__(self, hkl, experiment_id=0, UB=None):
        """
        Solve the prediction formula for the reflecting angle phi.

        If UB is given, override the contained crystal model. This is
        for use in refinement with time-varying crystal parameters
        """

        e = self._experiments[experiment_id]
        ray_predictor = ScanStaticRayPredictor(
            e.beam.get_s0(),
            e.goniometer.get_rotation_axis_datum(),
            e.goniometer.get_fixed_rotation(),
            e.goniometer.get_setting_rotation(),
            self._sequence_range,
        )

        UB_ = UB if UB else e.crystal.get_A()

        rays = ray_predictor(hkl, UB_)

        return rays


class ExperimentsPredictor:
    """
    Predict for relps based on the current states of models of the experimental
    geometry. This version manages multiple experiments, selecting the correct
    predictor in each case.
    """

    def __init__(self, experiments):
        """Construct by linking to instances of experimental model classes"""

        self._experiments = experiments

    def __call__(self, reflections):
        """Predict for all reflections at the current model geometry"""

        for iexp, e in enumerate(self._experiments):
            # select the reflections for this experiment only
            sel = reflections["id"] == iexp
            refs = reflections.select(sel)

            self._predict_one_experiment(e, refs)
            refs = self._post_predict_one_experiment(e, refs)

            # write predictions back to overall reflections
            reflections.set_selected(sel, refs)

        reflections = self._post_prediction(reflections)

        return reflections

    def _predict_one_experiment(self, experiment, reflections):
        raise NotImplementedError()

    def _post_predict_one_experiment(self, experiment, reflections):
        return reflections

    def _post_prediction(self, reflections):
        """Perform tasks on the whole reflection list after prediction before
        returning."""

        return reflections


class ScansExperimentsPredictor(ExperimentsPredictor):
    def _predict_one_experiment(self, experiment, reflections):
        # scan-varying
        if "ub_matrix" in reflections:
            predictor = sv(experiment)
            UB = reflections["ub_matrix"]
            s0 = reflections["s0_vector"]
            dmat = reflections["d_matrix"]
            Smat = reflections["S_matrix"]
            predictor.for_reflection_table(reflections, UB, s0, dmat, Smat)
        # scan static
        else:
            predictor = sc(experiment)
            UB = experiment.crystal.get_A()
            predictor.for_reflection_table(reflections, UB)

    def _post_prediction(self, reflections):
        if "xyzobs.mm.value" in reflections:
            reflections = self._match_full_turns(reflections)

        return reflections

    def _match_full_turns(self, reflections):
        """Modify the calculated phi values so that they match the full rotation
        from zero taken from the the observations, rather than being modulo 2*pi."""

        x_obs, y_obs, phi_obs = reflections["xyzobs.mm.value"].parts()
        x_calc, y_calc, phi_calc = reflections["xyzcal.mm"].parts()
        resid = phi_calc - (flex.fmod_positive(phi_obs, TWO_PI))
        # ensure this is the smaller of two possibilities
        resid = flex.fmod_positive((resid + pi), TWO_PI) - pi
        phi_calc = phi_obs + resid
        reflections["xyzcal.mm"] = flex.vec3_double(x_calc, y_calc, phi_calc)

        # Update xyzcal.px with the correct z_px values in keeping with above
        for iexp, e in enumerate(self._experiments):
            sel = reflections["id"] == iexp
            x_px, y_px, z_px = reflections["xyzcal.px"].select(sel).parts()
            scan = e.scan
            if scan is not None:
                z_px = scan.get_array_index_from_angle(phi_calc.select(sel), deg=False)
            else:
                # must be a still image, z centroid not meaningful
                z_px = phi_calc.select(sel)
            xyzcal_px = flex.vec3_double(x_px, y_px, z_px)
            reflections["xyzcal.px"].set_selected(sel, xyzcal_px)

        return reflections


class StillsExperimentsPredictor(ExperimentsPredictor):
    spherical_relp_model = False

    def _predict_one_experiment(self, experiment, reflections):
        predictor = st(experiment, spherical_relp=self.spherical_relp_model)
        UB = experiment.crystal.get_A()
        predictor.for_reflection_table(reflections, UB)


class LaueExperimentsPredictor(ExperimentsPredictor):
    def _predict_one_experiment(self, experiment, reflections):
        min_s0_idx = min(
            range(len(reflections["wavelength"])),
            key=reflections["wavelength"].__getitem__,
        )

        if "s0" not in reflections:
            unit_s0 = experiment.beam.get_unit_s0()
            wl = reflections["wavelength"][min_s0_idx]
            min_s0 = (unit_s0[0] / wl, unit_s0[1] / wl, unit_s0[2] / wl)
        else:
            min_s0 = reflections["s0"][min_s0_idx]

        dmin = experiment.detector.get_max_resolution(min_s0)
        predictor = LaueReflectionPredictor(experiment, dmin)
        UB = experiment.crystal.get_A()
        predictor.for_reflection_table(reflections, UB)


class TOFExperimentsPredictor(LaueExperimentsPredictor):
    def _post_predict_one_experiment(self, experiment, reflections):
        # Add ToF to xyzcal.mm
        wavelength_cal = reflections["wavelength_cal"]
        distance = experiment.beam.get_sample_to_source_distance() * 10**-3
        distance = distance + (reflections["s1"].norms() * 10**-3)
        tof_cal = tof_helpers.tof_from_wavelength(distance, wavelength_cal)  # (s)
        x, y, z = reflections["xyzcal.mm"].parts()
        tof_cal = tof_cal * 1e6  # (usec)
        reflections["xyzcal.mm"] = flex.vec3_double(x, y, tof_cal)

        # Add frame to xyzcal.px
        expt_tof = experiment.scan.get_property("time_of_flight")  # (usec)
        frames = list(range(len(expt_tof)))
        tof_to_frame = tof_helpers.tof_to_frame_interpolator(expt_tof, frames)
        tof_cal.set_selected(tof_cal < min(expt_tof), min(expt_tof))
        tof_cal.set_selected(tof_cal > max(expt_tof), max(expt_tof))
        reflection_frames = flex.double(tof_to_frame(tof_cal))
        px, py, pz = reflections["xyzcal.px"].parts()
        reflections["xyzcal.px"] = flex.vec3_double(px, py, reflection_frames)

        return reflections


class ExperimentsPredictorFactory:
    @staticmethod
    def from_experiments(experiments, force_stills=False, spherical_relp=False):
        assert experiments.all_same_type(), (
            "Cannot create ExperimentsPredictor for a mixture of experiments with different types"
        )

        if experiments.all_tof():
            return TOFExperimentsPredictor(experiments)
        elif experiments.all_laue():
            return LaueExperimentsPredictor(experiments)

        if not force_stills:
            for exp in experiments:
                if exp.goniometer is None:
                    force_stills = True
                    break

        if force_stills:
            predictor = StillsExperimentsPredictor(experiments)
            predictor.spherical_relp_model = spherical_relp

        else:
            predictor = ScansExperimentsPredictor(experiments)

        return predictor

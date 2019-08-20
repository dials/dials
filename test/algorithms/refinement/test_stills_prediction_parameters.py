from __future__ import absolute_import, division, print_function

import pytest

#### Python and general cctbx imports
from cctbx.sgtbx import space_group, space_group_symbols
from libtbx.phil import parse
from scitbx.math import five_number_summary

#### dials imports
from dials.array_family import flex

##### Import model builder
from dials.test.algorithms.refinement.setup_geometry import Extract

##### Imports for reflection prediction
from dials.algorithms.spot_prediction import IndexGenerator
from dxtbx.model.experiment_list import ExperimentList, Experiment
from dials.algorithms.refinement.prediction.managed_predictors import ScansRayPredictor

#### Import model parameterisations
from dials.algorithms.refinement.parameterisation.prediction_parameters_stills import (
    StillsPredictionParameterisation,
)
from dials.algorithms.refinement.parameterisation.prediction_parameters_stills import (
    SphericalRelpStillsPredictionParameterisation,
)

from dials.algorithms.refinement.parameterisation.detector_parameters import (
    DetectorParameterisationSinglePanel,
)
from dials.algorithms.refinement.parameterisation.beam_parameters import (
    BeamParameterisation,
)
from dials.algorithms.refinement.parameterisation.crystal_parameters import (
    CrystalOrientationParameterisation,
    CrystalUnitCellParameterisation,
)


class _Test(object):
    def __init__(self):
        self.create_models()
        self.generate_reflections()

    def create_models(self):
        # build models, with a larger crystal than default in order to get plenty of
        # reflections on the 'still' image
        overrides = """
    geometry.parameters.crystal.a.length.range=40 50;
    geometry.parameters.crystal.b.length.range=40 50;
    geometry.parameters.crystal.c.length.range=40 50;
    geometry.parameters.random_seed = 42"""

        master_phil = parse(
            """
        include scope dials.test.algorithms.refinement.geometry_phil
        """,
            process_includes=True,
        )

        models = Extract(master_phil, overrides)

        # keep track of the models
        self.detector = models.detector
        self.gonio = models.goniometer
        self.crystal = models.crystal
        self.beam = models.beam

        # Create a stills ExperimentList
        self.stills_experiments = ExperimentList()
        self.stills_experiments.append(
            Experiment(
                beam=self.beam,
                detector=self.detector,
                crystal=self.crystal,
                imageset=None,
            )
        )

        # keep track of the parameterisation of the models
        self.det_param = DetectorParameterisationSinglePanel(self.detector)
        self.s0_param = BeamParameterisation(self.beam, self.gonio)
        self.xlo_param = CrystalOrientationParameterisation(self.crystal)
        self.xluc_param = CrystalUnitCellParameterisation(self.crystal)

    def generate_reflections(self):
        # Build a mock scan for a 3 degree sweep
        from dxtbx.model import ScanFactory

        sf = ScanFactory()
        self.scan = sf.make_scan(
            image_range=(1, 1),
            exposure_times=0.1,
            oscillation=(0, 3.0),
            epochs=list(range(1)),
            deg=True,
        )
        sweep_range = self.scan.get_oscillation_range(deg=False)

        # Create a scans ExperimentList, only for generating reflections
        experiments = ExperimentList()
        experiments.append(
            Experiment(
                beam=self.beam,
                detector=self.detector,
                goniometer=self.gonio,
                scan=self.scan,
                crystal=self.crystal,
                imageset=None,
            )
        )

        # Create a ScansRayPredictor
        ray_predictor = ScansRayPredictor(experiments, sweep_range)

        # Generate rays - only to work out which hkls are predicted
        resolution = 2.0
        index_generator = IndexGenerator(
            self.crystal.get_unit_cell(),
            space_group(space_group_symbols(1).hall()).type(),
            resolution,
        )
        indices = index_generator.to_array()
        rays = ray_predictor(indices)

        # Make a standard reflection_table and copy in the ray data
        self.reflections = flex.reflection_table.empty_standard(len(rays))
        self.reflections.update(rays)

    def get_fd_gradients(self, pred_param, ref_predictor):

        # get finite difference gradients
        p_vals = pred_param.get_param_vals()
        deltas = [1.0e-7] * len(p_vals)

        fd_grads = []
        p_names = pred_param.get_param_names()
        for i in range(len(deltas)):

            # save parameter value
            val = p_vals[i]

            # calc reverse state
            p_vals[i] -= deltas[i] / 2.0
            pred_param.set_param_vals(p_vals)

            ref_predictor(self.reflections)

            x, y, _ = self.reflections["xyzcal.mm"].deep_copy().parts()
            delpsi = self.reflections["delpsical.rad"].deep_copy()
            rev_state = flex.vec3_double(x, y, delpsi)

            # calc forward state
            p_vals[i] += deltas[i]
            pred_param.set_param_vals(p_vals)

            ref_predictor(self.reflections)

            x, y, _ = self.reflections["xyzcal.mm"].deep_copy().parts()
            delpsi = self.reflections["delpsical.rad"].deep_copy()
            fwd_state = flex.vec3_double(x, y, delpsi)

            # reset parameter to saved value
            p_vals[i] = val

            # finite difference
            fd = fwd_state - rev_state
            x_grads, y_grads, delpsi_grads = fd.parts()
            x_grads /= deltas[i]
            y_grads /= deltas[i]
            delpsi_grads /= deltas[i]

            fd_grads.append(
                {
                    "name": p_names[i],
                    "dX_dp": x_grads,
                    "dY_dp": y_grads,
                    "dDeltaPsi_dp": delpsi_grads,
                }
            )

        # return to the initial state
        pred_param.set_param_vals(p_vals)

        return fd_grads


@pytest.fixture(scope="session")
def tc():
    return _Test()


def test_stills_pred_param(tc):
    print("Testing derivatives for StillsPredictionParameterisation")
    print("========================================================")

    # Build a prediction parameterisation for the stills experiment
    pred_param = StillsPredictionParameterisation(
        tc.stills_experiments,
        detector_parameterisations=[tc.det_param],
        beam_parameterisations=[tc.s0_param],
        xl_orientation_parameterisations=[tc.xlo_param],
        xl_unit_cell_parameterisations=[tc.xluc_param],
    )

    # Predict the reflections in place. Must do this ahead of calculating
    # the analytical gradients so quantities like s1 are correct
    from dials.algorithms.refinement.prediction.managed_predictors import (
        StillsExperimentsPredictor,
    )

    ref_predictor = StillsExperimentsPredictor(tc.stills_experiments)
    ref_predictor(tc.reflections)

    # get analytical gradients
    an_grads = pred_param.get_gradients(tc.reflections)

    fd_grads = tc.get_fd_gradients(pred_param, ref_predictor)

    for i, (an_grad, fd_grad) in enumerate(zip(an_grads, fd_grads)):

        # compare FD with analytical calculations
        print("\nParameter {0}: {1}".format(i, fd_grad["name"]))

        for name in ["dX_dp", "dY_dp", "dDeltaPsi_dp"]:
            print(name)
            a = fd_grad[name]
            b = an_grad[name]

            abs_error = a - b

            fns = five_number_summary(abs_error)
            print(
                ("  summary of absolute errors: %9.6f %9.6f %9.6f " + "%9.6f %9.6f")
                % fns
            )
            assert flex.max(flex.abs(abs_error)) < 0.0003
            # largest absolute error found to be about 0.00025 for dY/dp of
            # Crystal0g_param_3. Reject outlying absolute errors and test again.
            iqr = fns[3] - fns[1]

            # skip further stats on errors with an iqr of near zero, e.g. dDeltaPsi_dp
            # for detector parameters, which are all equal to zero
            if iqr < 1.0e-10:
                continue

            sel1 = abs_error < fns[3] + 1.5 * iqr
            sel2 = abs_error > fns[1] - 1.5 * iqr
            sel = sel1 & sel2
            tst = flex.max_index(flex.abs(abs_error.select(sel)))
            tst_val = abs_error.select(sel)[tst]
            n_outliers = sel.count(False)
            print(
                (
                    "  {0} outliers rejected, leaving greatest "
                    + "absolute error: {1:9.6f}"
                ).format(n_outliers, tst_val)
            )
            # largest absolute error now 0.000086 for dX/dp of Beam0Mu2
            assert abs(tst_val) < 0.00009

            # Completely skip parameters with FD gradients all zero (e.g. gradients of
            # DeltaPsi for detector parameters)
            sel1 = flex.abs(a) < 1.0e-10
            if sel1.all_eq(True):
                continue

            # otherwise calculate normalised errors, by dividing absolute errors by
            # the IQR (more stable than relative error calculation)
            norm_error = abs_error / iqr
            fns = five_number_summary(norm_error)
            print(
                ("  summary of normalised errors: %9.6f %9.6f %9.6f " + "%9.6f %9.6f")
                % fns
            )
            # largest normalised error found to be about 25.7 for dY/dp of
            # Crystal0g_param_3.
            try:
                assert flex.max(flex.abs(norm_error)) < 30
            except AssertionError as e:
                e.args += (
                    "extreme normalised error value: {}".format(
                        flex.max(flex.abs(norm_error))
                    ),
                )
                raise e

            # Reject outlying normalised errors and test again
            iqr = fns[3] - fns[1]
            if iqr > 0.0:
                sel1 = norm_error < fns[3] + 1.5 * iqr
                sel2 = norm_error > fns[1] - 1.5 * iqr
                sel = sel1 & sel2
                tst = flex.max_index(flex.abs(norm_error.select(sel)))
                tst_val = norm_error.select(sel)[tst]
                n_outliers = sel.count(False)

                # most outliers found for for dY/dp of Crystal0g_param_3 (which had
                # largest errors, so no surprise there).
                try:
                    assert n_outliers < 250
                except AssertionError as e:
                    e.args += ("too many outliers rejected: {}".format(n_outliers),)
                    raise e

                print(
                    (
                        "  {0} outliers rejected, leaving greatest "
                        + "normalised error: {1:9.6f}"
                    ).format(n_outliers, tst_val)
                )
                # largest normalied error now about -4. for dX/dp of Detector0Tau1
                assert abs(tst_val) < 6, "should be < 6, not %s" % tst_val


# In comparison with FD approximations, the worst gradients by far are dX/dp
# and dY/dp for parameter Crystal0g_param_3. Is this to do with the geometry
# of the test case?


def test_spherical_relp_stills_pred_param(tc):
    print("Testing derivatives for SphericalRelpStillsPredictionParameterisation")
    print("=====================================================================")

    # Build a prediction parameterisation for the stills experiment
    pred_param = SphericalRelpStillsPredictionParameterisation(
        tc.stills_experiments,
        detector_parameterisations=[tc.det_param],
        beam_parameterisations=[tc.s0_param],
        xl_orientation_parameterisations=[tc.xlo_param],
        xl_unit_cell_parameterisations=[tc.xluc_param],
    )

    # Predict the reflections in place. Must do this ahead of calculating
    # the analytical gradients so quantities like s1 are correct
    from dials.algorithms.refinement.prediction.managed_predictors import (
        ExperimentsPredictorFactory,
    )

    ref_predictor = ExperimentsPredictorFactory.from_experiments(
        tc.stills_experiments, force_stills=True, spherical_relp=True
    )
    ref_predictor(tc.reflections)

    # get analytical gradients
    an_grads = pred_param.get_gradients(tc.reflections)

    fd_grads = tc.get_fd_gradients(pred_param, ref_predictor)

    # compare FD with analytical calculations
    for i, (an_grad, fd_grad) in enumerate(zip(an_grads, fd_grads)):
        print("\nParameter {0}: {1}".format(i, fd_grad["name"]))
        for name in ["dX_dp", "dY_dp", "dDeltaPsi_dp"]:
            print(name)
            for a, b in zip(an_grad[name], fd_grad[name]):
                if name == "dDeltaPsi_dp":
                    # DeltaPsi errors are much worse than X, Y errors!
                    # FIXME, look into this further
                    assert a == pytest.approx(b, abs=5e-3)
                else:
                    assert a == pytest.approx(b, abs=5e-6)
            print("OK")

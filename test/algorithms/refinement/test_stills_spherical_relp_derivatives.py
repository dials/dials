"""
Test derivatives typed up in dials_regression/doc/notes/prediction/stills_prediction_nave3.pdf
"""

from __future__ import absolute_import, division, print_function

import pytest

from cctbx.sgtbx import space_group, space_group_symbols
from libtbx.phil import parse
from scitbx import matrix
from dials.array_family import flex
from dials.test.algorithms.refinement.setup_geometry import Extract
from dials.algorithms.spot_prediction import IndexGenerator
from dxtbx.model.experiment_list import ExperimentList, Experiment
from dials.algorithms.refinement.prediction.managed_predictors import ScansRayPredictor
from dials.algorithms.refinement.parameterisation.prediction_parameters_stills import (
    StillsPredictionParameterisation,
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

# Create a reflection predictor specific for this test
from dials.algorithms.spot_prediction import StillsReflectionPredictor


class Predictor(object):
    def __init__(self, experiments):

        self._experiment = experiments[0]
        self.update()

    def update(self):
        """Build predictor objects for the current geometry of each Experiment"""

        self._predictor = StillsReflectionPredictor(
            self._experiment, spherical_relp=True
        )
        self._UB = matrix.sqr(self._experiment.crystal.get_U()) * matrix.sqr(
            self._experiment.crystal.get_B()
        )

    def predict(self, reflections):
        """
        Predict for all reflections
        """

        self._predictor.for_reflection_table(reflections, self._UB)
        return reflections


# Simplied equivalent of a PredictionParameterisation class for this test
class AnalyticalGradients(object):
    """A class to implement the analytical gradient calculation in the document"""

    def __init__(
        self,
        experiments,
        detector_parameterisation,
        beam_parameterisation,
        xl_orientation_parameterisation,
        xl_unit_cell_parameterisation,
    ):

        # References to the underlying models from the first experiment
        self.experiment = experiments[0]
        self.beam = self.experiment.beam
        self.crystal = self.experiment.crystal

        # Keep references to all parameterisations
        self.detector_parameterisation = detector_parameterisation
        self.beam_parameterisation = beam_parameterisation
        self.xl_orientation_parameterisation = xl_orientation_parameterisation
        self.xl_unit_cell_parameterisation = xl_unit_cell_parameterisation

        # Keep quantities of interest from the model
        self.U = matrix.sqr(self.crystal.get_U())
        self.B = matrix.sqr(self.crystal.get_B())
        self.UB = self.U * self.B
        self.s0 = matrix.col(self.beam.get_s0())
        self.s0len = self.s0.length()
        self.s0len2 = self.s0.length_sq()
        self.us0 = self.s0.normalize()

    def get_beam_gradients(self, reflections):

        ds0_dbeam_p = self.beam_parameterisation.get_ds_dp()
        p_names = self.beam_parameterisation.get_param_names()

        n = len(reflections)
        U = flex.mat3_double(n, self.U)
        B = flex.mat3_double(n, self.B)
        UB = U * B

        # q is the reciprocal lattice vector, in the lab frame
        h = reflections["miller_index"].as_vec3_double()
        q = UB * h
        qlen2 = q.dot(q)

        q_s0 = q + self.s0
        s1 = reflections["s1"]
        ss = qlen2 + 2 * q.dot(self.s0) + self.s0len2
        assert (ss > 0.0).all_eq(True)
        s = flex.sqrt(ss)
        sss = s * ss
        inv_s = 1.0 / s
        inv_sss = 1.0 / sss

        # check equation 10
        tmp = self.s0len * (q_s0) / s
        for a, b in zip(s1, tmp):
            assert a == pytest.approx(b, abs=1e-7)

        ds1_dp = {}

        # loop through the parameters
        for name, der in zip(p_names, ds0_dbeam_p):

            # term1
            term1 = self.us0.dot(der) * q_s0 + self.s0len * (der)
            term1 = term1 * inv_s

            # term2
            term2 = self.s0len * q_s0 * q_s0.dot(der)
            term2 = term2 * inv_sss

            name = "Beam1" + name  # XXXX Hack to get matching keys
            ds1_dp[name] = {"ds1": (term1 - term2)}

        return ds1_dp

    def get_crystal_orientation_gradients(self, reflections):

        # get derivatives of the U matrix wrt the parameters
        dU_dxlo_p = self.xl_orientation_parameterisation.get_ds_dp()
        p_names = self.xl_orientation_parameterisation.get_param_names()

        n = len(reflections)
        U = flex.mat3_double(n, self.U)
        B = flex.mat3_double(n, self.B)
        UB = U * B

        # q is the reciprocal lattice vector, in the lab frame
        h = reflections["miller_index"].as_vec3_double()
        q = UB * h
        qlen2 = q.dot(q)

        q_s0 = q + self.s0
        ss = qlen2 + 2 * q.dot(self.s0) + self.s0len2
        assert (ss > 0.0).all_eq(True)
        s = flex.sqrt(ss)
        sss = s * ss
        inv_s = 1.0 / s
        inv_sss = 1.0 / sss

        ds1_dp = {}

        # loop through the parameters
        for name, der in zip(p_names, dU_dxlo_p):

            # calculate the derivative of q for this parameter
            dq = flex.mat3_double(n, der.elems) * B * h

            # term1
            term1 = self.s0len * dq
            term1 = term1 * inv_s

            # term2
            term2 = self.s0len * q_s0 * q_s0.dot(dq)
            term2 = term2 * inv_sss

            name = "Crystal1" + name  # XXXX Hack to get matching keys
            ds1_dp[name] = {"ds1": (term1 - term2)}

        return ds1_dp

    def get_crystal_unit_cell_gradients(self, reflections):

        # get derivatives of the B matrix wrt the parameters
        dB_dxluc_p = self.xl_unit_cell_parameterisation.get_ds_dp()
        p_names = self.xl_unit_cell_parameterisation.get_param_names()

        n = len(reflections)
        U = flex.mat3_double(n, self.U)
        B = flex.mat3_double(n, self.B)
        UB = U * B

        # q is the reciprocal lattice vector, in the lab frame
        h = reflections["miller_index"].as_vec3_double()
        q = UB * h
        qlen2 = q.dot(q)

        q_s0 = q + self.s0
        ss = qlen2 + 2 * q.dot(self.s0) + self.s0len2
        assert (ss > 0.0).all_eq(True)
        s = flex.sqrt(ss)
        sss = s * ss
        inv_s = 1.0 / s
        inv_sss = 1.0 / sss

        ds1_dp = {}

        # loop through the parameters
        for name, der in zip(p_names, dB_dxluc_p):

            # calculate the derivative of q for this parameter
            dq = U * flex.mat3_double(n, der.elems) * h

            # term1
            term1 = self.s0len * dq
            term1 = term1 * inv_s

            # term2
            term2 = self.s0len * q_s0 * q_s0.dot(dq)
            term2 = term2 * inv_sss

            name = "Crystal1" + name  # XXXX Hack to get matching keys
            ds1_dp[name] = {"ds1": (term1 - term2)}

        return ds1_dp


def test():
    # Build models, with a larger crystal than default in order to get plenty of
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

    mydetector = models.detector
    mygonio = models.goniometer
    mycrystal = models.crystal
    mybeam = models.beam

    # Build a mock scan for a 3 degree sweep
    from dxtbx.model import ScanFactory

    sf = ScanFactory()
    myscan = sf.make_scan(
        image_range=(1, 1),
        exposure_times=0.1,
        oscillation=(0, 3.0),
        epochs=list(range(1)),
        deg=True,
    )
    sweep_range = myscan.get_oscillation_range(deg=False)

    # Create parameterisations of these models
    det_param = DetectorParameterisationSinglePanel(mydetector)
    s0_param = BeamParameterisation(mybeam, mygonio)
    xlo_param = CrystalOrientationParameterisation(mycrystal)
    xluc_param = CrystalUnitCellParameterisation(mycrystal)

    # Create a scans ExperimentList, only for generating reflections
    experiments = ExperimentList()
    experiments.append(
        Experiment(
            beam=mybeam,
            detector=mydetector,
            goniometer=mygonio,
            scan=myscan,
            crystal=mycrystal,
            imageset=None,
        )
    )

    # Create a stills ExperimentList
    stills_experiments = ExperimentList()
    stills_experiments.append(
        Experiment(beam=mybeam, detector=mydetector, crystal=mycrystal, imageset=None)
    )

    # Generate rays - only to work out which hkls are predicted
    ray_predictor = ScansRayPredictor(experiments, sweep_range)
    resolution = 2.0
    index_generator = IndexGenerator(
        mycrystal.get_unit_cell(),
        space_group(space_group_symbols(1).hall()).type(),
        resolution,
    )
    indices = index_generator.to_array()
    rays = ray_predictor(indices)

    # Make a standard reflection_table and copy in the ray data
    reflections = flex.reflection_table.empty_standard(len(rays))
    reflections.update(rays)

    # Build a standard prediction parameterisation for the stills experiment to do
    # FD calculation (not used for its analytical gradients)
    pred_param = StillsPredictionParameterisation(
        stills_experiments,
        detector_parameterisations=[det_param],
        beam_parameterisations=[s0_param],
        xl_orientation_parameterisations=[xlo_param],
        xl_unit_cell_parameterisations=[xluc_param],
    )

    # Make a managed SphericalRelpStillsReflectionPredictor reflection predictor
    # for the first (only) experiment
    ref_predictor = Predictor(stills_experiments)

    # Predict these reflections in place. Must do this ahead of calculating
    # the analytical gradients so quantities like s1 are correct
    ref_predictor.update()
    ref_predictor.predict(reflections)

    # calculate analytical gradients
    ag = AnalyticalGradients(
        stills_experiments,
        detector_parameterisation=det_param,
        beam_parameterisation=s0_param,
        xl_orientation_parameterisation=xlo_param,
        xl_unit_cell_parameterisation=xluc_param,
    )
    an_grads = ag.get_beam_gradients(reflections)
    an_grads.update(ag.get_crystal_orientation_gradients(reflections))
    an_grads.update(ag.get_crystal_unit_cell_gradients(reflections))

    # get finite difference gradients
    p_vals = pred_param.get_param_vals()
    deltas = [1.0e-7] * len(p_vals)

    fd_grads = []
    p_names = pred_param.get_param_names()
    for i, delta in enumerate(deltas):

        # save parameter value
        val = p_vals[i]

        # calc reverse state
        p_vals[i] -= delta / 2.0
        pred_param.set_param_vals(p_vals)

        ref_predictor.update()
        ref_predictor.predict(reflections)

        x, y, _ = reflections["xyzcal.mm"].deep_copy().parts()
        s1 = reflections["s1"].deep_copy()
        rev_state = s1

        # calc forward state
        p_vals[i] += delta
        pred_param.set_param_vals(p_vals)

        ref_predictor.update()
        ref_predictor.predict(reflections)

        x, y, _ = reflections["xyzcal.mm"].deep_copy().parts()
        s1 = reflections["s1"].deep_copy()
        fwd_state = s1

        # reset parameter to saved value
        p_vals[i] = val

        # finite difference - currently for s1 only
        fd = fwd_state - rev_state
        inv_delta = 1.0 / delta
        s1_grads = fd * inv_delta

        # store gradients
        fd_grads.append({"name": p_names[i], "ds1": s1_grads})

    # return to the initial state
    pred_param.set_param_vals(p_vals)

    for i, fd_grad in enumerate(fd_grads):

        ## compare FD with analytical calculations
        print("\n\nParameter {0}: {1}".format(i, fd_grad["name"]))

        print("d[s1]/dp for the first reflection")
        print("finite diff", fd_grad["ds1"][0])
        try:
            an_grad = an_grads[fd_grad["name"]]
        except KeyError:
            continue

        print("checking analytical vs finite difference gradients for s1")
        for a, b in zip(fd_grad["ds1"], an_grad["ds1"]):
            assert a == pytest.approx(b, abs=1e-7)

from __future__ import annotations

import random
import sys
import time
from math import pi

import numpy as np
import pytest

import dxtbx.flumpy as flumpy
from dxtbx.model.experiment_list import Experiment, ExperimentList
from scitbx.array_family import flex

from dials.algorithms.refinement.parameterisation.scan_varying_beam_parameters import (
    ScanVaryingBeamParameterisation,
)
from dials.algorithms.refinement.parameterisation.scan_varying_crystal_parameters import (
    ScanVaryingCrystalOrientationParameterisation,
    ScanVaryingCrystalUnitCellParameterisation,
)
from dials.algorithms.refinement.parameterisation.scan_varying_detector_parameters import (
    ScanVaryingDetectorParameterisationSinglePanel,
)
from dials.algorithms.refinement.parameterisation.scan_varying_goniometer_parameters import (
    ScanVaryingGoniometerParameterisation,
)
from dials.algorithms.refinement.parameterisation.scan_varying_prediction_parameters import (
    ScanVaryingPredictionParameterisation,
    SparseFlex,
)
from dials.algorithms.refinement.prediction.managed_predictors import (
    ScansExperimentsPredictor,
    ScansRayPredictor,
)
from dials_refinement_helpers_ext import (
    build_reconstitute_derivatives_mat3,
    build_reconstitute_derivatives_vec3,
    intersection_i_seqs_unsorted,
)


class _Test:
    def create_models(self, cmdline_overrides=None):
        from dxtbx.model import ScanFactory
        from libtbx.phil import parse

        from . import geometry_phil
        from .setup_geometry import Extract

        if cmdline_overrides is None:
            cmdline_overrides = []
        overrides = """geometry.parameters.crystal.a.length.range = 10 50
geometry.parameters.crystal.b.length.range = 10 50
geometry.parameters.crystal.c.length.range = 10 50"""

        master_phil = parse(geometry_phil)

        # Extract models
        models = Extract(master_phil, overrides, cmdline_args=cmdline_overrides)
        self.detector = models.detector
        self.goniometer = models.goniometer
        self.crystal = models.crystal
        self.beam = models.beam

        # Make a scan of 1-20 * 0.5 deg images
        sf = ScanFactory()
        self.scan = sf.make_scan((1, 20), 0.5, (0, 0.5), list(range(20)))

        # Generate an ExperimentList
        self.experiments = ExperimentList()
        self.experiments.append(
            Experiment(
                beam=self.beam,
                detector=self.detector,
                goniometer=self.goniometer,
                scan=self.scan,
                crystal=self.crystal,
                imageset=None,
            )
        )

        # Create a reflection predictor for the experiments
        self.ref_predictor = ScansExperimentsPredictor(self.experiments)

        # Create scan-varying parameterisations of these models, with 3 samples
        self.det_param = ScanVaryingDetectorParameterisationSinglePanel(
            self.detector, self.scan.get_array_range(), 3
        )
        self.s0_param = ScanVaryingBeamParameterisation(
            self.beam, self.scan.get_array_range(), 3, self.goniometer
        )
        self.xlo_param = ScanVaryingCrystalOrientationParameterisation(
            self.crystal, self.scan.get_array_range(), 3
        )
        self.xluc_param = ScanVaryingCrystalUnitCellParameterisation(
            self.crystal, self.scan.get_array_range(), 3
        )
        self.gon_param = ScanVaryingGoniometerParameterisation(
            self.goniometer, self.scan.get_array_range(), 3, self.beam
        )

    def generate_reflections(self):
        from cctbx.sgtbx import space_group, space_group_symbols

        from dials.algorithms.spot_prediction import IndexGenerator, ray_intersection

        sequence_range = self.scan.get_oscillation_range(deg=False)
        resolution = 2.0
        index_generator = IndexGenerator(
            self.crystal.get_unit_cell(),
            space_group(space_group_symbols(1).hall()).type(),
            resolution,
        )
        indices = index_generator.to_array()

        # Predict rays within the sequence range
        ray_predictor = ScansRayPredictor(self.experiments, sequence_range)
        obs_refs = ray_predictor(indices)

        # Take only those rays that intersect the detector
        intersects = ray_intersection(self.detector, obs_refs)
        obs_refs = obs_refs.select(intersects)

        # Re-predict using the Experiments predictor for all these reflections. The
        # result is the same, but we gain also the flags and xyzcal.px columns
        obs_refs["id"] = flex.int(len(obs_refs), 0)
        obs_refs = self.ref_predictor(obs_refs)

        # Set 'observed' centroids from the predicted ones
        obs_refs["xyzobs.mm.value"] = obs_refs["xyzcal.mm"]

        # Invent some variances for the centroid positions of the simulated data
        im_width = 0.1 * pi / 180.0
        px_size = self.detector[0].get_pixel_size()
        var_x = flex.double(len(obs_refs), (px_size[0] / 2.0) ** 2)
        var_y = flex.double(len(obs_refs), (px_size[1] / 2.0) ** 2)
        var_phi = flex.double(len(obs_refs), (im_width / 2.0) ** 2)
        obs_refs["xyzobs.mm.variance"] = flex.vec3_double(var_x, var_y, var_phi)

        # set the flex random seed to an 'uninteresting' number
        flex.set_random_seed(12407)

        # take 10 random reflections for speed
        reflections = obs_refs.select(flex.random_selection(len(obs_refs), 10))

        # use a BlockCalculator to calculate the blocks per image
        from dials.algorithms.refinement.reflection_manager import BlockCalculator

        block_calculator = BlockCalculator(self.experiments, reflections)
        reflections = block_calculator.per_image()

        return reflections


def test(cmdline_overrides=[]):
    tc = _Test()
    tc.create_models(cmdline_overrides)
    reflections = tc.generate_reflections()

    # use a ReflectionManager to exclude reflections too close to the spindle,
    # plus set the frame numbers
    from dials.algorithms.refinement.reflection_manager import ReflectionManager

    refman = ReflectionManager(reflections, tc.experiments, outlier_detector=None)
    refman.finalise()

    # create prediction parameterisation of the requested type
    pred_param = ScanVaryingPredictionParameterisation(
        tc.experiments,
        [tc.det_param],
        [tc.s0_param],
        [tc.xlo_param],
        [tc.xluc_param],
        [tc.gon_param],
    )

    # keep only those reflections that pass inclusion criteria and have predictions
    reflections = refman.get_matches()

    # get analytical gradients
    pred_param.compose(reflections)
    an_grads = pred_param.get_gradients(reflections)

    # get finite difference gradients
    p_vals = pred_param.get_param_vals()
    p_names = pred_param.get_param_names()
    deltas = [1.0e-7] * len(p_vals)

    for i, delta in enumerate(deltas):
        val = p_vals[i]

        p_vals[i] -= delta / 2.0
        pred_param.set_param_vals(p_vals)
        pred_param.compose(reflections)

        tc.ref_predictor(reflections)

        rev_state = reflections["xyzcal.mm"].deep_copy()

        p_vals[i] += delta
        pred_param.set_param_vals(p_vals)
        pred_param.compose(reflections)

        tc.ref_predictor(reflections)

        fwd_state = reflections["xyzcal.mm"].deep_copy()
        p_vals[i] = val

        fd = fwd_state - rev_state
        x_grads, y_grads, phi_grads = fd.parts()
        x_grads /= delta
        y_grads /= delta
        phi_grads /= delta

        try:
            for (a, b) in zip(x_grads, an_grads[i]["dX_dp"]):
                assert a == pytest.approx(b, abs=1e-5)
            for (a, b) in zip(y_grads, an_grads[i]["dY_dp"]):
                assert a == pytest.approx(b, abs=1e-5)
            for (a, b) in zip(phi_grads, an_grads[i]["dphi_dp"]):
                assert a == pytest.approx(b, abs=1e-5)
        except AssertionError:
            print(f"Failure for {p_names[i]}")
            raise

    # return to the initial state
    pred_param.set_param_vals(p_vals)
    pred_param.compose(reflections)


def test_SparseFlex_scalars():

    size = 100

    # Make a dense double array with 50% explicit zeroes
    arr = flex.random_double(size)
    indices = flex.random_selection(size, int(size / 2))
    elements = arr.select(indices)
    arr *= 0.0
    arr.set_selected(indices, elements)

    # Make the equivalent SparseFlex
    sf_arr = SparseFlex(size, elements, indices)

    # Test multiplication of SparseFlex[double] by scalar
    sf2 = sf_arr * 2.0  # __mul__
    for a, b in zip(sf2.as_dense_vector(), arr * 2.0):
        assert a == b

    sf2 = 2.0 * sf_arr  # __rmul__
    for a, b in zip(sf2.as_dense_vector(), 2.0 * arr):
        assert a == b

    # Test division of SparseFlex[double] by scalar
    sf_half = sf_arr / 2.0
    for a, b in zip(sf_half.as_dense_vector(), arr / 2.0):
        assert a == b

    # Test addition of two SparseFlex[double]s
    sf2 = sf_arr + sf_arr
    for a, b in zip(sf2.as_dense_vector(), arr + arr):
        assert a == b

    # Test subtraction between two SparseFlex[double]s
    sf3 = sf2 - sf_arr
    for a, b in zip(sf3.as_dense_vector(), arr):
        assert a == b


def test_SparseFlex_matrix_and_vector_arithmetic():

    size = 100

    # Make a dense vec3 array with 50% explicit zeroes
    vec = flex.vec3_double(
        flex.random_double(size), flex.random_double(size), flex.random_double(size)
    )
    indices = flex.random_selection(size, int(size / 2))
    elements = vec.select(indices)
    vec *= 0.0
    vec.set_selected(indices, elements)

    # Make the equivalent SparseFlex[vec3]
    sf_vec = SparseFlex(size, elements, indices)

    # Make a dense mat3 array with the same 50% explicit zeroes
    mat = flex.mat3_double(
        (flex.random_double_r3_rotation_matrix() for i in range(size))
    )
    elements = mat.select(indices)
    mat *= 0.0
    mat.set_selected(indices, elements)

    # Make the equivalent SparseFlex[mat3]
    sf_mat = SparseFlex(size, elements, indices)

    # Test multiplication of SparseFlex[vec3] by scalar
    sf2 = sf_vec * 2.0  # __mul__
    for a, b in zip(sf2.as_dense_vector(), vec * 2.0):
        assert a == b

    sf2 = 2.0 * sf_vec  # __rmul__
    for a, b in zip(sf2.as_dense_vector(), vec * 2.0):
        assert a == b

    # Test multiplication of SparseFlex[mat3] by scalar. Only __mul__, because
    # mat3_double does not have __rmul__ defined
    sf2 = sf_mat * 2.0
    for a, b in zip(sf2.as_dense_vector(), mat * 2.0):
        assert a == b

    # Test matrix multiplication: SparseFlex[mat3] * flex.vec3_double. Use a
    # new vector which does not have explicit zero elements
    vec2 = flex.vec3_double(
        flex.random_double(size), flex.random_double(size), flex.random_double(size)
    )
    sf_rot = sf_mat * vec2
    for a, b in zip(sf_rot.as_dense_vector(), mat * vec2):
        assert a == b

    # Test matrix multiplication: SparseFlex[mat3] * SparseFlex[vec3]
    sf_rot = sf_mat * sf_vec
    for a, b in zip(sf_rot.as_dense_vector(), mat * vec):
        assert a == b

    # Test matrix multiplication: SparseFlex[mat3] * flex.mat3_double. Use a
    # new matrix which does not have explicit zero elements
    mat2 = flex.mat3_double(
        (flex.random_double_r3_rotation_matrix() for i in range(size))
    )
    sf_mat2 = sf_mat * mat2
    for a, b in zip(sf_mat2.as_dense_vector(), mat * mat2):
        assert a == b

    # Test matrix multiplication SparseFlex[mat3] * SparseFlex[mat3]
    sf_mat2 = sf_mat * sf_mat
    for a, b in zip(sf_mat2.as_dense_vector(), mat * mat):
        assert a == b

    # Test addition of two SparseFlex[vec3]s
    sf2 = sf_vec + sf_vec
    for a, b in zip(sf2.as_dense_vector(), vec + vec):
        assert a == b

    # Test subtraction between two SparseFlex[vec3]s
    sf3 = sf2 - sf_vec
    for a, b in zip(sf3.as_dense_vector(), vec):
        assert a == b


def test_SparseFlex_vec3_only_methods():

    size = 100

    # Make a dense vec3 array with 50% explicit zeroes
    vec = flex.vec3_double(
        flex.random_double(size), flex.random_double(size), flex.random_double(size)
    )
    indices = flex.random_selection(size, int(size / 2))
    elements = vec.select(indices)
    vec *= 0.0
    vec.set_selected(indices, elements)

    # Make the equivalent SparseFlex[vec3]
    sf_vec = SparseFlex(size, elements, indices)

    # Test dot product of SparseFlex[vec3] and flex.vec3_double. Use a
    # new vector which does not have explicit zero elements
    vec2 = flex.vec3_double(
        flex.random_double(size), flex.random_double(size), flex.random_double(size)
    )
    dotprod = sf_vec.dot(vec2)
    for a, b in zip(dotprod.as_dense_vector(), vec.dot(vec2)):
        assert a == b

    # Test dot product of SparseFlex[vec3] and SparseFlex[vec3]
    dotprod = sf_vec.dot(sf_vec)
    for a, b in zip(dotprod.as_dense_vector(), vec.dot(vec)):
        assert a == b

    # Test rotate_around_origin
    directions = flex.vec3_double(
        flex.random_double(size), flex.random_double(size), flex.random_double(size)
    )
    angles = flex.random_double(size) * pi
    sf_rot = sf_vec.rotate_around_origin(directions, angles)
    for a, b in zip(
        sf_rot.as_dense_vector(), vec.rotate_around_origin(directions, angles)
    ):
        assert a == b

    # Test parts
    sf_x, sf_y, sf_z = sf_vec.parts()
    x, y, z = vec.parts()
    for a, b in zip(sf_x.as_dense_vector(), x):
        assert a == b
    for a, b in zip(sf_y.as_dense_vector(), y):
        assert a == b
    for a, b in zip(sf_z.as_dense_vector(), z):
        assert a == b


def test_SparseFlex_select():

    size = 100

    # Make a dense double array with 50% explicit zeroes
    arr = flex.random_double(size)
    indices = flex.random_selection(size, int(size / 2))
    elements = arr.select(indices)
    arr *= 0.0
    arr.set_selected(indices, elements)

    # Make the equivalent SparseFlex
    sf_arr = SparseFlex(size, elements, indices)

    # Make a new random selection
    isel = flex.random_selection(size, int(size / 2))

    # Select
    sf2 = sf_arr.select(isel)
    arr2 = arr.select(isel)

    # Check
    for a, b in zip(sf2.as_dense_vector(), arr2):
        assert a == b


@pytest.mark.parametrize(
    "random_order", [False, pytest.param(True, marks=pytest.mark.xfail)]
)
def test_SparseFlex_select_intersection(random_order):
    """SparseFlex.select requires the intersection of two size_t arrays, A and
    B. We view this as the selection of elements in A according to their
    presence in B. We have three approaches for this: a pure Python method,
    a C++ equivalent and a Numpy method. The Numpy method is fast, but fails
    when the input arrays are not sorted. We cannot assume this to be the case
    in practice, so have to reject that method for use in SparseFlex.select.

    With changes to the C++ version, this is now about 8 times faster than the
    pure Python version and 4 times faster than the NumPy version."""

    # Two random selections
    size = 1000
    sel1 = flex.random_selection(size, int(size / 2))
    if random_order:
        perm = flex.random_permutation(int(size / 2))
        sel1 = sel1.select(perm)
    sel2 = flex.random_selection(size, int(size / 2))
    if random_order:
        perm = flex.random_permutation(int(size / 2))
        sel2 = sel2.select(perm)

    # Pure Python version
    start = time.perf_counter_ns()
    index_a = flex.size_t(0)
    index_b = flex.size_t(0)
    lookup = {}
    for i_a, val in enumerate(sel1):
        lookup[val] = i_a
    for i_b, val in enumerate(sel2):
        i_a = lookup.get(val)
        if i_a is not None:
            index_a.append(i_a)
            index_b.append(i_b)
    end = time.perf_counter_ns()
    wc_time_py = end - start

    # C++ version
    start = time.perf_counter_ns()
    index_a_cpp, index_b_cpp = intersection_i_seqs_unsorted(sel1, sel2)
    end = time.perf_counter_ns()
    wc_time_cpp = end - start

    # Numpy version
    start = time.perf_counter_ns()
    _, index_a_np, index_b_np = np.intersect1d(
        flumpy.to_numpy(sel1),
        flumpy.to_numpy(sel2),
        assume_unique=True,
        return_indices=True,
    )
    index_a_np = flumpy.from_numpy(index_a_np)
    index_b_np = flumpy.from_numpy(index_b_np)
    end = time.perf_counter_ns()
    wc_time_numpy = end - start

    print(f"Timing info Python: {wc_time_py} C++ {wc_time_cpp} Numpy: {wc_time_numpy}")

    # Check results are equal: Python to C++
    for (
        a,
        b,
    ) in zip(index_a_cpp, index_a):
        assert a == b
    for (
        a,
        b,
    ) in zip(index_b_cpp, index_b):
        assert a == b

    # Check results are equal: Python to Numpy
    for (
        a,
        b,
    ) in zip(index_a_np, index_a):
        assert a == b
    for (
        a,
        b,
    ) in zip(index_b_np, index_b):
        assert a == b


def test_intersection_i_seqs_speed():

    exec_times = []
    for i in range(100):

        size = 10000

        sel1 = flex.random_selection(size, int(size / 2))
        sel2 = flex.random_selection(size, int(size / 2))

        # C++ version
        start = time.perf_counter_ns()
        index_a_cpp, index_b_cpp = intersection_i_seqs_unsorted(sel1, sel2)
        end = time.perf_counter_ns()
        wc_time_cpp = end - start

        # Python version
        start = time.perf_counter_ns()
        index_a = flex.size_t(0)
        index_b = flex.size_t(0)
        lookup = {}
        for i_a, val in enumerate(sel1):
            lookup[val] = i_a
        for i_b, val in enumerate(sel2):
            i_a = lookup.get(val)
            if i_a is not None:
                index_a.append(i_a)
                index_b.append(i_b)
        end = time.perf_counter_ns()
        wc_time_py = end - start

        for (
            a,
            b,
        ) in zip(index_a_cpp, index_a):
            assert a == b
        for (
            a,
            b,
        ) in zip(index_b_cpp, index_b):
            assert a == b

        exec_times.append((wc_time_cpp, wc_time_py))

    tot_cpp, tot_py = zip(*exec_times)
    print(f"Total time in C++ function: {sum(tot_cpp)}")
    print(f"Total time in Python version: {sum(tot_py)}")

    # We use the C++ version because it is faster than the Python version.
    # Let's ensure that's the case across platforms
    assert sum(tot_cpp) / sum(tot_py) < 1.0


@pytest.mark.parametrize("element_type", ["vec3", "mat3"])
def test_ReconstituteDerivatives(element_type):

    if element_type == "vec3":
        n = 3
        build = build_reconstitute_derivatives_vec3
        flex_type = flex.vec3_double
    else:
        n = 9
        build = build_reconstitute_derivatives_mat3
        flex_type = flex.mat3_double

    # Create 100 random (element, indices) pairs, with up to 100 elements in
    # the indices
    pair_data = []
    total_size = 0
    for i in range(100):
        element = [random.random() for j in range(n)]
        indices = flex.random_size_t(random.randint(0, 100))
        pair_data.append((element, indices))
        total_size += indices.size()

    rec_der = build(total_size)

    reference_data = flex_type(0)
    reference_indices = flex.size_t(0)
    for element, indices in pair_data:
        rec_der.add_data(element, indices)
        n = indices.size()
        values = flex_type(n, element)
        reference_data.extend(values)
        reference_indices.extend(indices)

    # Compare the reference values from extending flex arrays with the
    # reconstituted values
    test_data = rec_der.get_data()
    test_indices = rec_der.get_indices()

    assert test_data.as_double().all_eq(reference_data.as_double())
    assert test_indices.all_eq(reference_indices)


if __name__ == "__main__":
    cmdline_overrides = sys.argv[1:]
    test(cmdline_overrides)

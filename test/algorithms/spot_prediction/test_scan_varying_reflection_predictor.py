from __future__ import absolute_import, division, print_function

import copy
import os
import pytest


class Data(object):
    def __init__(self, dials_regression):
        path = os.path.join(
            dials_regression,
            "prediction_test_data",
            "experiments_scan_varying_crystal.json",
        )

        from dxtbx.model.experiment_list import ExperimentListFactory

        self.experiments = ExperimentListFactory.from_json_file(path)
        assert len(self.experiments) == 1
        assert (
            self.experiments[0].crystal.num_scan_points
            == self.experiments[0].scan.get_num_images() + 1
        )

    def run(self):
        self.tst_regression()
        # self.tst_with_hkl()
        # self.tst_with_hkl_and_panel()
        # self.tst_with_hkl_and_panel_list()
        self.tst_for_reflection_table()
        self.tst_vs_static()

    def _predict_new(self, hkl=None, frame=None, panel=None):
        from dials.algorithms.spot_prediction import ScanVaryingReflectionPredictor
        from dials.array_family import flex

        predict = ScanVaryingReflectionPredictor(self.experiments[0])
        # if hkl is None:
        A = [
            self.experiments[0].crystal.get_A_at_scan_point(i)
            for i in range(self.experiments[0].crystal.num_scan_points)
        ]
        result = predict.for_ub(flex.mat3_double(A))
        # else:
        # if panel is None:
        # result = predict(hkl, frame)
        # else:
        # result = predict(hkl, frame, panel)
        return result


@pytest.fixture(scope="session")
def data(dials_regression):
    return Data(dials_regression)


def test_regression(data):
    r_new = data._predict_new()
    assert len(r_new) == 1934
    assert r_new[0]["miller_index"] == (-8, -30, -23)
    assert r_new[0]["xyzcal.px"] == pytest.approx(
        (75.33831543451907, 2327.55737978813, 0.2548567552525226), abs=1e-6
    )

    # def tst_with_hkl(self):
    # from dials.algorithms.spot_prediction import ReekeIndexGenerator
    # from dials.array_family import flex
    # from scitbx import matrix
    # import scitbx

    # m2 = self.experiments[0].goniometer.get_rotation_axis()
    # s0 = self.experiments[0].beam.get_s0()
    # dmin = self.experiments[0].detector.get_max_resolution(s0)
    # margin = 1
    # scan = self.experiments[0].scan
    # crystal = self.experiments[0].crystal
    # frame_0 = scan.get_array_range()[0]
    # step = 1

    # all_indices = flex.miller_index()
    # all_frames = flex.int()
    # for frame in range(*scan.get_array_range()):

    # phi_beg = scan.get_angle_from_array_index(frame, deg = False)
    # phi_end = scan.get_angle_from_array_index(frame + step, deg = False)
    # r_beg = matrix.sqr(scitbx.math.r3_rotation_axis_and_angle_as_matrix(
    # axis = m2, angle = phi_beg, deg = False))
    # r_end = matrix.sqr(scitbx.math.r3_rotation_axis_and_angle_as_matrix(
    # axis = m2, angle = phi_end, deg = False))

    # A1 = r_beg * crystal.get_A_at_scan_point(frame - frame_0)

    # A2 = r_end * crystal.get_A_at_scan_point(frame - frame_0 + step)

    # indices = ReekeIndexGenerator(A1, A2, m2, s0, dmin, margin)
    # indices = indices.to_array()
    # all_indices.extend(indices)
    # all_frames.extend(flex.int(len(indices), frame))

    # r_old = self._predict_new()
    # r_new = self._predict_new(all_indices, all_frames)
    # assert(len(r_old) == len(r_new))
    # print 'OK'

    # def tst_with_hkl_and_panel(self):
    # from dials.algorithms.spot_prediction import ReekeIndexGenerator
    # from dials.array_family import flex
    # from scitbx import matrix
    # import scitbx

    # m2 = self.experiments[0].goniometer.get_rotation_axis()
    # s0 = self.experiments[0].beam.get_s0()
    # dmin = self.experiments[0].detector.get_max_resolution(s0)
    # margin = 1
    # scan = self.experiments[0].scan
    # crystal = self.experiments[0].crystal
    # frame_0 = scan.get_array_range()[0]
    # step = 1

    # all_indices = flex.miller_index()
    # all_frames = flex.int()
    # for frame in range(*scan.get_array_range()):

    # phi_beg = scan.get_angle_from_array_index(frame, deg = False)
    # phi_end = scan.get_angle_from_array_index(frame + step, deg = False)
    # r_beg = matrix.sqr(scitbx.math.r3_rotation_axis_and_angle_as_matrix(
    # axis = m2, angle = phi_beg, deg = False))
    # r_end = matrix.sqr(scitbx.math.r3_rotation_axis_and_angle_as_matrix(
    # axis = m2, angle = phi_end, deg = False))

    # A1 = r_beg * crystal.get_A_at_scan_point(frame - frame_0)

    # A2 = r_end * crystal.get_A_at_scan_point(frame - frame_0 + step)

    # indices = ReekeIndexGenerator(A1, A2, m2, s0, dmin, margin)
    # indices = indices.to_array()
    # all_indices.extend(indices)
    # all_frames.extend(flex.int(len(indices), frame))

    # r_old = self._predict_new()
    # try:
    # r_new = self._predict_new(all_indices, all_frames, 1)
    # assert(False)
    # except Exception:
    # pass

    # r_new = self._predict_new(all_indices, all_frames, 0)
    # assert(len(r_old) < len(r_new))
    # print 'OK'

    # def tst_with_hkl_and_panel_list(self):

    # from dials.algorithms.spot_prediction import ReekeIndexGenerator
    # from dials.array_family import flex
    # from scitbx import matrix
    # import scitbx

    # m2 = self.experiments[0].goniometer.get_rotation_axis()
    # s0 = self.experiments[0].beam.get_s0()
    # dmin = self.experiments[0].detector.get_max_resolution(s0)
    # margin = 1
    # scan = self.experiments[0].scan
    # crystal = self.experiments[0].crystal
    # frame_0 = scan.get_array_range()[0]
    # step = 1

    # all_indices = flex.miller_index()
    # all_frames = flex.int()
    # for frame in range(*scan.get_array_range()):

    # phi_beg = scan.get_angle_from_array_index(frame, deg = False)
    # phi_end = scan.get_angle_from_array_index(frame + step, deg = False)
    # r_beg = matrix.sqr(scitbx.math.r3_rotation_axis_and_angle_as_matrix(
    # axis = m2, angle = phi_beg, deg = False))
    # r_end = matrix.sqr(scitbx.math.r3_rotation_axis_and_angle_as_matrix(
    # axis = m2, angle = phi_end, deg = False))

    # A1 = r_beg * crystal.get_A_at_scan_point(frame - frame_0)

    # A2 = r_end * crystal.get_A_at_scan_point(frame - frame_0 + step)

    # indices = ReekeIndexGenerator(A1, A2, m2, s0, dmin, margin)
    # indices = indices.to_array()
    # all_indices.extend(indices)
    # all_frames.extend(flex.int(len(indices), frame))

    # r_old = self._predict_new()
    # try:
    # r_new = self._predict_new(all_indices, all_frames,
    # flex.size_t(len(all_indices), 1))
    # assert(False)
    # except Exception:
    # pass

    # r_new = self._predict_new(all_indices, all_frames,
    # flex.size_t(len(all_indices), 0))
    # assert(len(r_old) < len(r_new))
    # print 'OK'


from dials.test.algorithms.spot_prediction.test_scan_static_reflection_predictor import (  # noqa: F401, used as test fixture
    data as static_test,
)


def test_scan_varying_results_are_close_to_static_prediction_when_model_is_static(
    static_test  # noqa: F811, not a redefinition
):
    """Test that various modes of scan-varying prediction produce results
    close to static prediction when the supplied models are indeed static"""

    # Get static predictor results

    scan = static_test.experiments[0].scan
    crystal = static_test.experiments[0].crystal
    beam = static_test.experiments[0].beam
    goniometer = static_test.experiments[0].goniometer
    n_scan_points = scan.get_num_images() + 1
    static_preds = static_test.predict_new()
    static_preds.sort("miller_index")

    # Set up scan-varying predictor
    from dials.algorithms.spot_prediction import ScanVaryingReflectionPredictor
    from dials.array_family import flex

    predict = ScanVaryingReflectionPredictor(static_test.experiments[0])

    def compare(refs1, refs2):
        assert len(refs1) == len(refs2)
        for r1, r2 in zip(refs1, refs2):
            assert r1["miller_index"] == r2["miller_index"]
            # differences less than one hundredth of a pixel/image
            for e1, e2 in zip(r1["xyzcal.px"], r2["xyzcal.px"]):
                assert e1 == pytest.approx(e2, abs=0.01)

    # Prediction for UB matrix expressed as array of static UB
    A = [crystal.get_A() for i in range(n_scan_points)]
    result1 = predict.for_ub(flex.mat3_double(A))
    result1.sort("miller_index")
    compare(static_preds, result1)

    # Prediction for UB matrix, s0 vectors and goniometer setting rotation
    # matrices expressed as arrays of static model states
    s0 = [beam.get_s0() for i in range(n_scan_points)]
    S = [goniometer.get_setting_rotation() for i in range(n_scan_points)]
    result2 = predict.for_varying_models(
        flex.mat3_double(A), flex.vec3_double(s0), flex.mat3_double(S)
    )
    result2.sort("miller_index")
    compare(static_preds, result2)

    # First frame only, start and end UB
    _, _, z = static_preds["xyzcal.px"].parts()
    static_preds_frame0 = static_preds.select((z >= 0) & (z < 1))
    A = crystal.get_A()
    result3 = predict.for_ub_on_single_image(0, A, A)
    result3.sort("miller_index")
    compare(static_preds_frame0, result3)

    # First frame only, start and end UB, s0 and S
    s0 = beam.get_s0()
    S = goniometer.get_setting_rotation()
    result4 = predict.for_varying_models_on_single_image(0, A, A, s0, s0, S, S)
    result4.sort("miller_index")
    compare(static_preds_frame0, result4)


def test_for_reflection_table(data):
    from dials.algorithms.spot_prediction import (
        ScanVaryingReflectionPredictor,
        ScanStaticReflectionPredictor,
    )
    from dials.array_family import flex

    predict = ScanStaticReflectionPredictor(data.experiments[0])
    preds = predict.for_ub(data.experiments[0].crystal.get_A())

    preds["ub_matrix"] = flex.mat3_double(
        len(preds), data.experiments[0].crystal.get_A()
    )
    preds["s0"] = flex.vec3_double(len(preds), data.experiments[0].beam.get_s0())
    preds["d_matrix"] = flex.mat3_double(len(preds))
    preds["S_matrix"] = flex.mat3_double(
        len(preds), data.experiments[0].goniometer.get_setting_rotation()
    )
    for ipanel, panel in enumerate(data.experiments[0].detector):
        sel = preds["panel"] == ipanel
        D = panel.get_d_matrix()
        preds["d_matrix"].set_selected(sel, D)
    predict = ScanVaryingReflectionPredictor(data.experiments[0])
    old_preds = copy.deepcopy(preds)
    predict.for_reflection_table(
        preds, preds["ub_matrix"], preds["s0"], preds["d_matrix"], preds["S_matrix"]
    )

    # Because UB, s0, d and S values are the same for all reflections, the new
    # reflections should be approx equal to those produced by the scan static
    # predictor
    old_x, old_y, old_z = old_preds["xyzcal.px"].parts()
    new_x, new_y, new_z = preds["xyzcal.px"].parts()
    assert old_x.all_approx_equal(new_x)
    assert old_y.all_approx_equal(new_y)
    assert old_z.all_approx_equal(new_z)

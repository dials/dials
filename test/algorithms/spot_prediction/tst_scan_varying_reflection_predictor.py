from __future__ import absolute_import, division

class Test(object):

  def __init__(self):
    import libtbx.load_env
    try:
      dials_regression = libtbx.env.dist_path('dials_regression')
    except KeyError, e:
      print 'SKIP: dials_regression not configured'
      exit(0)

    import os
    path = os.path.join(
      dials_regression,
      'prediction_test_data',
      'experiments_scan_varying_crystal.json')

    from dxtbx.model.experiment_list import ExperimentListFactory
    self.experiments = ExperimentListFactory.from_json_file(path)
    assert(len(self.experiments) == 1)
    assert(self.experiments[0].crystal.num_scan_points ==
           self.experiments[0].scan.get_num_images() + 1)

  def run(self):
    self.tst_regression()
    #self.tst_with_hkl()
    #self.tst_with_hkl_and_panel()
    #self.tst_with_hkl_and_panel_list()
    self.tst_for_reflection_table()

  def predict_new(self, hkl=None, frame=None, panel=None):
    from dials.algorithms.spot_prediction import ScanVaryingReflectionPredictor
    from time import time
    from dials.array_family import flex
    st = time()
    predict = ScanVaryingReflectionPredictor(self.experiments[0])
    #if hkl is None:
    A = [self.experiments[0].crystal.get_A_at_scan_point(i) for i in
           range(self.experiments[0].crystal.num_scan_points)]
    result = predict.for_ub(flex.mat3_double(A))
    #else:
      #if panel is None:
        #result = predict(hkl, frame)
      #else:
        #result = predict(hkl, frame, panel)

    #print "New Time: ", time() - st
    return result

  def tst_regression(self):
    from libtbx.test_utils import approx_equal
    r_new = self.predict_new()
    assert len(r_new) == 1934
    assert r_new[0]['miller_index'] == (-8, -30, -23)
    assert approx_equal(r_new[0]['xyzcal.px'],
      (75.33831543451907, 2327.55737978813, 0.2548567552525226))
    print "OK"
    return

  #def tst_with_hkl(self):
    #from dials.algorithms.spot_prediction import ReekeIndexGenerator
    #from dials.array_family import flex
    #from scitbx import matrix
    #import scitbx

    #m2 = self.experiments[0].goniometer.get_rotation_axis()
    #s0 = self.experiments[0].beam.get_s0()
    #dmin = self.experiments[0].detector.get_max_resolution(s0)
    #margin = 1
    #scan = self.experiments[0].scan
    #crystal = self.experiments[0].crystal
    #frame_0 = scan.get_array_range()[0]
    #step = 1

    #all_indices = flex.miller_index()
    #all_frames = flex.int()
    #for frame in range(*scan.get_array_range()):

      #phi_beg = scan.get_angle_from_array_index(frame, deg = False)
      #phi_end = scan.get_angle_from_array_index(frame + step, deg = False)
      #r_beg = matrix.sqr(scitbx.math.r3_rotation_axis_and_angle_as_matrix(
          #axis = m2, angle = phi_beg, deg = False))
      #r_end = matrix.sqr(scitbx.math.r3_rotation_axis_and_angle_as_matrix(
          #axis = m2, angle = phi_end, deg = False))

      #A1 = r_beg * crystal.get_A_at_scan_point(frame - frame_0)

      #A2 = r_end * crystal.get_A_at_scan_point(frame - frame_0 + step)

      #indices = ReekeIndexGenerator(A1, A2, m2, s0, dmin, margin)
      #indices = indices.to_array()
      #all_indices.extend(indices)
      #all_frames.extend(flex.int(len(indices), frame))

    #r_old = self.predict_new()
    #r_new = self.predict_new(all_indices, all_frames)
    #assert(len(r_old) == len(r_new))
    #print 'OK'

  #def tst_with_hkl_and_panel(self):
    #from dials.algorithms.spot_prediction import ReekeIndexGenerator
    #from dials.array_family import flex
    #from scitbx import matrix
    #import scitbx

    #m2 = self.experiments[0].goniometer.get_rotation_axis()
    #s0 = self.experiments[0].beam.get_s0()
    #dmin = self.experiments[0].detector.get_max_resolution(s0)
    #margin = 1
    #scan = self.experiments[0].scan
    #crystal = self.experiments[0].crystal
    #frame_0 = scan.get_array_range()[0]
    #step = 1

    #all_indices = flex.miller_index()
    #all_frames = flex.int()
    #for frame in range(*scan.get_array_range()):

      #phi_beg = scan.get_angle_from_array_index(frame, deg = False)
      #phi_end = scan.get_angle_from_array_index(frame + step, deg = False)
      #r_beg = matrix.sqr(scitbx.math.r3_rotation_axis_and_angle_as_matrix(
          #axis = m2, angle = phi_beg, deg = False))
      #r_end = matrix.sqr(scitbx.math.r3_rotation_axis_and_angle_as_matrix(
          #axis = m2, angle = phi_end, deg = False))

      #A1 = r_beg * crystal.get_A_at_scan_point(frame - frame_0)

      #A2 = r_end * crystal.get_A_at_scan_point(frame - frame_0 + step)

      #indices = ReekeIndexGenerator(A1, A2, m2, s0, dmin, margin)
      #indices = indices.to_array()
      #all_indices.extend(indices)
      #all_frames.extend(flex.int(len(indices), frame))

    #r_old = self.predict_new()
    #try:
      #r_new = self.predict_new(all_indices, all_frames, 1)
      #assert(False)
    #except Exception:
      #pass

    #r_new = self.predict_new(all_indices, all_frames, 0)
    #assert(len(r_old) < len(r_new))
    #print 'OK'

  #def tst_with_hkl_and_panel_list(self):

    #from dials.algorithms.spot_prediction import ReekeIndexGenerator
    #from dials.array_family import flex
    #from scitbx import matrix
    #import scitbx

    #m2 = self.experiments[0].goniometer.get_rotation_axis()
    #s0 = self.experiments[0].beam.get_s0()
    #dmin = self.experiments[0].detector.get_max_resolution(s0)
    #margin = 1
    #scan = self.experiments[0].scan
    #crystal = self.experiments[0].crystal
    #frame_0 = scan.get_array_range()[0]
    #step = 1

    #all_indices = flex.miller_index()
    #all_frames = flex.int()
    #for frame in range(*scan.get_array_range()):

      #phi_beg = scan.get_angle_from_array_index(frame, deg = False)
      #phi_end = scan.get_angle_from_array_index(frame + step, deg = False)
      #r_beg = matrix.sqr(scitbx.math.r3_rotation_axis_and_angle_as_matrix(
          #axis = m2, angle = phi_beg, deg = False))
      #r_end = matrix.sqr(scitbx.math.r3_rotation_axis_and_angle_as_matrix(
          #axis = m2, angle = phi_end, deg = False))

      #A1 = r_beg * crystal.get_A_at_scan_point(frame - frame_0)

      #A2 = r_end * crystal.get_A_at_scan_point(frame - frame_0 + step)

      #indices = ReekeIndexGenerator(A1, A2, m2, s0, dmin, margin)
      #indices = indices.to_array()
      #all_indices.extend(indices)
      #all_frames.extend(flex.int(len(indices), frame))

    #r_old = self.predict_new()
    #try:
      #r_new = self.predict_new(all_indices, all_frames,
                               #flex.size_t(len(all_indices), 1))
      #assert(False)
    #except Exception:
      #pass

    #r_new = self.predict_new(all_indices, all_frames,
                             #flex.size_t(len(all_indices), 0))
    #assert(len(r_old) < len(r_new))
    #print 'OK'

  def tst_for_reflection_table(self):
    from libtbx.test_utils import approx_equal
    from dials.algorithms.spot_prediction import \
      ScanVaryingReflectionPredictor, ScanStaticReflectionPredictor
    from dials.array_family import flex

    predict = ScanStaticReflectionPredictor(self.experiments[0])
    preds = predict.for_ub(self.experiments[0].crystal.get_A())

    preds['ub_matrix'] = flex.mat3_double(len(preds), self.experiments[0].crystal.get_A())
    preds['s0'] = flex.vec3_double(len(preds), self.experiments[0].beam.get_s0())
    preds['d_matrix'] = flex.mat3_double(len(preds))
    for ipanel, panel in enumerate(self.experiments[0].detector):
      sel = preds['panel'] == ipanel
      D = panel.get_d_matrix()
      preds['d_matrix'].set_selected(sel, D)
    predict = ScanVaryingReflectionPredictor(self.experiments[0])
    from copy import deepcopy
    old_preds = deepcopy(preds)
    predict.for_reflection_table(preds,
                                 preds['ub_matrix'],
                                 preds['s0'],
                                 preds['d_matrix'])

    # Because UB, s0 and d values are the same for all reflections, the new
    # reflections should be approx equal to those produced by the scan static
    # predictor
    old_x, old_y, old_z = old_preds['xyzcal.px'].parts()
    new_x, new_y, new_z = preds['xyzcal.px'].parts()
    assert old_x.all_approx_equal(new_x)
    assert old_y.all_approx_equal(new_y)
    assert old_z.all_approx_equal(new_z)

    print "OK"

    return

if __name__ == '__main__':

  from dials.test import cd_auto
  with cd_auto(__file__):
    test = Test()
    test.run()

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
    from dxtbx.model.experiment_list import ExperimentListFactory
    path = os.path.join(
      dials_regression,
      'centroid_test_data',
      'experiments.json')

    self.experiments = ExperimentListFactory.from_json_file(path)
    assert(len(self.experiments) == 1)
    self.experiments[0].imageset.set_beam(self.experiments[0].beam)
    self.experiments[0].imageset.set_detector(self.experiments[0].detector)
    self.experiments[0].imageset.set_goniometer(self.experiments[0].goniometer)
    self.experiments[0].imageset.set_scan(self.experiments[0].scan)

    reflection_filename = os.path.join(
      dials_regression,
      'prediction_test_data',
      'expected_reflections.pickle')

    from dials.array_family import flex
    self.reflections = flex.reflection_table.from_pickle(reflection_filename)

  def run(self):
    self.tst_number_of_predictions()
    self.tst_vs_old()
    self.tst_with_reflection_table()
    #self.tst_with_hkl()
    #self.tst_with_hkl_and_panel()
    #self.tst_with_hkl_and_panel_list()

  def predict_new(self, experiment=None, hkl=None, panel=None):
    from dials.algorithms.spot_prediction import ScanStaticReflectionPredictor
    if experiment is None: experiment=self.experiments[0]
    predict = ScanStaticReflectionPredictor(experiment)
    if hkl is None:
      return predict.for_ub(experiment.crystal.get_A())
    else:
      if panel is None:
        return predict(hkl)
      else:
        return predict(hkl, panel)

  def tst_number_of_predictions(self):
    prediction_count = len(self.predict_new())
    assert prediction_count == 1996, "Predicted %d spots != 1996" % prediction_count

    import copy
    modified_experiment = copy.deepcopy(self.experiments[0])
    modified_experiment.scan.set_oscillation((360, 0.2))
    prediction_count = len(self.predict_new(experiment=modified_experiment))
    assert prediction_count == 1996, "Predicted %d spots != 1996 when rotated by 360 degrees" % prediction_count

  def tst_vs_old(self):
    r_old = self.reflections
    r_new = self.predict_new()
    assert(len(r_old) == len(r_new))
    eps = 1e-7
    for r1, r2 in zip(r_old.rows(), r_new.rows()):
      assert(r1['miller_index'] == r2['miller_index'])
      assert(r1['panel'] == r2['panel'])
      assert(r1['entering'] == r2['entering'])
      assert(all(abs(a-b) < eps for a, b in zip(r1['s1'], r2['s1'])))
      assert(all(abs(a-b) < eps for a, b in zip(r1['xyzcal.px'], r2['xyzcal.px'])))
      assert(all(abs(a-b) < eps for a, b in zip(r1['xyzcal.mm'], r2['xyzcal.mm'])))
    print 'OK'

  def tst_with_reflection_table(self):
    from dials.algorithms.spot_prediction import ScanStaticReflectionPredictor
    from dials.array_family import flex
    r_old = self.reflections
    r_new = flex.reflection_table()
    r_new['miller_index'] = r_old['miller_index']
    r_new['panel'] = r_old['panel']
    r_new['entering'] = r_old['entering']
    predict = ScanStaticReflectionPredictor(self.experiments[0])
    predict.for_reflection_table(r_new, self.experiments[0].crystal.get_A())
    eps = 1e-7
    for r1, r2 in zip(r_old.rows(), r_new.rows()):
      assert(r1['miller_index'] == r2['miller_index'])
      assert(r1['panel'] == r2['panel'])
      assert(r1['entering'] == r2['entering'])
      assert(all(abs(a-b) < eps for a, b in zip(r1['s1'], r2['s1'])))
      assert(all(abs(a-b) < eps for a, b in zip(r1['xyzcal.px'], r2['xyzcal.px'])))
      assert(all(abs(a-b) < eps for a, b in zip(r1['xyzcal.mm'], r2['xyzcal.mm'])))
    print 'OK'

  # FIXME Fix the Spot prediction interface
  #def tst_with_hkl(self):
    #from dials.algorithms.spot_prediction import IndexGenerator

    #unit_cell = self.experiments[0].crystal.get_unit_cell()
    #space_group_type = self.experiments[0].crystal.get_space_group().type()
    #dmin = self.experiments[0].detector.get_max_resolution(
      #self.experiments[0].beam.get_s0())
    #indices = IndexGenerator(unit_cell, space_group_type, dmin)

    #r_old = self.predict_new()
    #r_new = self.predict_new(self.experiments[0].crystal.get_A(),indices.to_array())
    #assert(len(r_old) == len(r_new))
    #eps = 1e-7
    #for r1, r2 in zip(r_old.rows(), r_new.rows()):
      #assert(r1['miller_index'] == r2['miller_index'])
      #assert(r1['panel'] == r2['panel'])
      #assert(r1['entering'] == r2['entering'])
      #assert(all(abs(a-b) < eps for a, b in zip(r1['s1'], r2['s1'])))
      #assert(all(abs(a-b) < eps for a, b in zip(r1['xyzcal.px'], r2['xyzcal.px'])))
      #assert(all(abs(a-b) < eps for a, b in zip(r1['xyzcal.mm'], r2['xyzcal.mm'])))
    #print 'OK'

  #def tst_with_hkl_and_panel(self):
    #from dials.algorithms.spot_prediction import IndexGenerator

    #unit_cell = self.experiments[0].crystal.get_unit_cell()
    #space_group_type = self.experiments[0].crystal.get_space_group().type()
    #dmin = self.experiments[0].detector.get_max_resolution(
      #self.experiments[0].beam.get_s0())
    #indices = IndexGenerator(unit_cell, space_group_type, dmin)
    #indices = indices.to_array()
    #try:
      #r_new = self.predict_new(indices, 1)
      #assert(False)
    #except Exception:
      #pass

    #r_old = self.predict_new()
    #r_new = self.predict_new(indices,0)
    #assert(len(r_old) < len(r_new))
    #print 'OK'

  #def tst_with_hkl_and_panel_list(self):
    #from dials.algorithms.spot_prediction import IndexGenerator
    #from dials.array_family import flex
    #unit_cell = self.experiments[0].crystal.get_unit_cell()
    #space_group_type = self.experiments[0].crystal.get_space_group().type()
    #dmin = self.experiments[0].detector.get_max_resolution(
      #self.experiments[0].beam.get_s0())
    #indices = IndexGenerator(unit_cell, space_group_type, dmin)
    #indices = indices.to_array()

    #panels = flex.size_t(len(indices), 1)
    #try:
      #r_new = self.predict_new(indices, panels)
      #assert(False)
    #except Exception:
      #pass

    #panels = flex.size_t(len(indices), 0)
    #r_old = self.predict_new()
    #r_new = self.predict_new(indices, panels)
    #assert(len(r_old) < len(r_new))
    #print 'OK'
    #pass


if __name__ == '__main__':
  from dials.test import cd_auto
  with cd_auto(__file__):

    test = Test()
    test.run()

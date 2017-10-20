
from __future__ import division

def run():
  from dials.algorithms.spot_prediction import PixelToMillerIndex
  from dials.array_family import flex
  from dxtbx.model.experiment_list import ExperimentListFactory
  import libtbx.load_env
  from os.path import join
  try:
    dials_regression = libtbx.env.dist_path('dials_regression')
  except KeyError:
    print 'SKIP: dials_regression not configured'
    exit(0)


  filename = join(dials_regression, "centroid_test_data/experiments.json")

  experiments = ExperimentListFactory.from_json_file(filename)

  transform = PixelToMillerIndex(
    experiments[0].beam,
    experiments[0].detector,
    experiments[0].goniometer,
    experiments[0].scan,
    experiments[0].crystal)

  reflections = flex.reflection_table.from_predictions(experiments[0])

  for r in reflections:
    panel = r['panel']
    x, y, z = r['xyzcal.px']
    h0 = r['miller_index']
    h1 = transform.h(panel, x, y, z)
    try:
      assert abs(h0[0] - h1[0]) < 1e-7
      assert abs(h0[1] - h1[1]) < 1e-7
      assert abs(h0[2] - h1[2]) < 1e-7
    except Exception:
      print h0, h1
      raise

  print 'OK'

if __name__ == '__main__':

  run()

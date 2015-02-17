
from __future__ import division

def run_single(experiments1):
  from dials.util.nexus import dump, load

  try:
    run_single.count += 1
  except Exception:
    run_single.count = 0

  # Dump the file
  dump(experiments1, None, "hklout_%d.nxs" % run_single.count)

  # Load the file
  experiments2, reflections = load("hklout_%d.nxs" % run_single.count)
  assert(experiments2 is not None)
  assert(reflections is None)

  print 'OK'

def run():
  from dxtbx.model.experiment.experiment_list import ExperimentListFactory
  from os.path import join
  import libtbx.load_env
  try:
    dials_regression = libtbx.env.dist_path('dials_regression')
  except KeyError, e:
    print 'FAIL: dials_regression not configured'
    exit(0)
  path = join(dials_regression, "nexus_test_data", "shared_models")
  filename_list = [
    'single.json',
    'multiple_unrelated.json',
    'multi_crystal.json',
    'two_colour.json',
    'multiple_sweeps.json',
    'stills.json'
  ]
  for filename in filename_list:
    experiments = ExperimentListFactory.from_json_file(join(path, filename))
    print filename
    run_single(experiments)

if __name__ == '__main__':
  from dials.test import cd_auto
  with cd_auto(__file__):
    run()

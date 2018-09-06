from __future__ import absolute_import, division, print_function

import six.moves.cPickle as pickle
import os
import procrunner

def test_import_integrate_hkl(dials_regression, tmpdir):
  tmpdir.chdir()

  from dials.array_family import flex # import dependency

  result = procrunner.run_process([
      'dials.import_xds',
      'input.method=reflections',
      os.path.join(dials_regression, "centroid_test_data", 'INTEGRATE.HKL'),
      os.path.join(dials_regression, "centroid_test_data", "experiments.json"),
  ])
  assert result['exitcode'] == 0
  assert result['stderr'] == ''

  with open('integrate_hkl.pickle', 'rb') as fh:
    table = pickle.load(fh)

  assert 'miller_index' in table
  assert 'id' in table
  assert 'panel' in table
  assert 'xyzcal.px' in table
  assert 'xyzobs.px.value' in table
  assert 'intensity.cor.value' in table
  assert 'intensity.cor.variance' in table
  assert len(table) == 174911

def test_import_spot_xds(dials_regression, tmpdir):
  tmpdir.chdir()

  result = procrunner.run_process([
      'dials.import_xds',
      'input.method=reflections',
      os.path.join(dials_regression, "centroid_test_data", 'SPOT.XDS'),
  ])
  assert result['exitcode'] == 0
  assert result['stderr'] == ''

  with open('spot_xds.pickle', 'rb') as fh:
    table = pickle.load(fh)

  assert 'miller_index' in table
  assert 'id' in table
  assert 'panel' in table
  assert 'xyzobs.px.value' in table
  assert 'intensity.sum.value' in table
  assert len(table) == 742

def test_import_spot_xds_with_filtering(dials_regression, tmpdir):
  tmpdir.chdir()

  result = procrunner.run_process([
      'dials.import_xds',
      'input.method=reflections',
      os.path.join(dials_regression, "centroid_test_data", 'SPOT.XDS'),
      'remove_invalid=True',
  ])
  assert result['exitcode'] == 0
  assert result['stderr'] == ''

  with open('spot_xds.pickle', 'rb') as fh:
    table = pickle.load(fh)

  assert 'miller_index' in table
  assert 'id' in table
  assert 'panel' in table
  assert 'xyzobs.px.value' in table
  assert 'intensity.sum.value' in table
  assert len(table) == 664

def test_from_xds_files(dials_regression, tmpdir):
  tmpdir.chdir()

  # Import from the image files
  result = procrunner.run_process([
      'dials.import_xds',
      'input.method=experiment',
      'output.filename=import_experiments.json',
      os.path.join(dials_regression, "centroid_test_data"),
  ])
  assert result['exitcode'] == 0
  assert result['stderr'] == ''

  assert os.path.exists("import_experiments.json")

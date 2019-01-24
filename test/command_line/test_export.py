from __future__ import absolute_import, division, print_function

import json
import os
import procrunner
import pytest
from dxtbx.serialize.load import _decode_dict

# Tests used to check for h5py
# May need to add this again if lack of this check causes issues.

def test_nxs(dials_regression, tmpdir):
  tmpdir.chdir()

  # Call dials.export
  result = procrunner.run_process([
      'dials.export',
      'format=nxs',
      os.path.join(dials_regression, "centroid_test_data", "experiments.json"),
      os.path.join(dials_regression, "centroid_test_data", "integrated.pickle"),
  ])
  assert result['exitcode'] == 0
  assert result['stderr'] == ''
  assert os.path.exists("integrated.nxs")

def test_mtz(dials_regression, tmpdir):
  tmpdir.chdir()

  # Call dials.export
  result = procrunner.run_process([
      'dials.export',
      'format=mtz',
      os.path.join(dials_regression, "centroid_test_data", "experiments.json"),
      os.path.join(dials_regression, "centroid_test_data", "integrated.pickle"),
  ])
  assert result['exitcode'] == 0
  assert result['stderr'] == ''
  assert os.path.exists("integrated.mtz")

def test_mmcif(dials_regression, tmpdir):
  tmpdir.chdir()

  # Call dials.export after integration
  result = procrunner.run_process([
      'dials.export',
      'format=mmcif',
      os.path.join(dials_regression, "centroid_test_data", "experiments.json"),
      os.path.join(dials_regression, "centroid_test_data", "integrated.pickle"),
  ])
  assert result['exitcode'] == 0
  assert result['stderr'] == ''
  assert os.path.exists("integrated.cif")

  #TODO include similar test for exporting scaled data in mmcif format

def test_xds_ascii(dials_regression, tmpdir):
  tmpdir.chdir()

  # Call dials.export
  result = procrunner.run_process([
      'dials.export',
      'intensity=sum',
      'format=xds_ascii',
      os.path.join(dials_regression, "centroid_test_data", "experiments.json"),
      os.path.join(dials_regression, "centroid_test_data", "integrated.pickle"),
  ])
  assert result['exitcode'] == 0
  assert result['stderr'] == ''
  assert os.path.exists("DIALS.HKL")

  psi_values = {
    (-9, 7, -10):153.430361,
    (-5, 11, -26):175.559441,
    (-4, 23, 24):129.468070,
    (2, 10, 20):147.947274
    }

  for record in open('DIALS.HKL', 'r'):
    if record.startswith('!'):
      continue
    tokens = record.split()
    hkl = tuple(map(int, tokens[:3]))
    if not hkl in psi_values:
      continue
    psi = float(tokens[-1])
    assert psi == pytest.approx(psi_values[hkl], abs=0.1)

def test_sadabs(dials_regression, tmpdir):
  tmpdir.chdir()

  # Call dials.export
  result = procrunner.run_process([
    'dials.export',
    'intensity=sum',
    'mtz.partiality_threshold=0.99',
    'format=sadabs',
      os.path.join(dials_regression, "centroid_test_data", "experiments.json"),
      os.path.join(dials_regression, "centroid_test_data", "integrated.pickle"),
  ])
  assert result['exitcode'] == 0
  assert result['stderr'] == ''
  assert os.path.exists("integrated.sad")

  direction_cosines = {
    (-9, 7, -10):(0.51253, -0.72107, 0.84696, -0.68476, -0.14130, -0.10561),
    (-5, 11, -26):(0.51310, -0.62895, 0.84711, -0.59223, -0.13830, -0.50366),
    (-4, 23, 24):(0.51308, -0.60578, 0.84711, -0.31416, -0.13840, 0.73099),
    (2, 10, 20):(0.51239, -0.46605, 0.84693, -0.61521, -0.14204, 0.63586)
    }

  with open('integrated.sad', 'r') as fh:
    for record in fh:
      record = record.replace('-', ' -')
      tokens = record.split()
      hkl = tuple(map(int, tokens[:3]))
      cosines = tuple(map(float, tokens[6:12]))
      if not hkl in direction_cosines:
        continue
      assert cosines == pytest.approx(direction_cosines[hkl], abs=0.001)

def test_json(dials_regression, tmpdir):
  tmpdir.chdir()

  # Call dials.export
  result = procrunner.run_process([
      'dials.export',
      'format=json',
      os.path.join(dials_regression, "centroid_test_data", "datablock.json"),
      os.path.join(dials_regression, "centroid_test_data", "strong.pickle"),
  ])
  assert result['exitcode'] == 0
  assert result['stderr'] == ''
  assert os.path.exists('rlp.json')

  from dxtbx.model.experiment_list import ExperimentListFactory
  with open('rlp.json', 'rb') as f:
    d = json.load(f, object_hook=_decode_dict)
  assert d.keys() == ['imageset_id', 'experiments', 'rlp', 'experiment_id'], d.keys()
  assert d['rlp'][:3] == [0.123454, 0.57687, 0.186465], d['rlp'][:3]
  assert d['imageset_id'][0] == 0
  assert d['experiment_id'][0] == 0
  experiments = ExperimentListFactory.from_dict(d['experiments'])
  imgset = experiments.imagesets()
  assert len(imgset) == 1

def test_json_shortened(dials_regression, tmpdir):
  tmpdir.chdir()

  # Call dials.export
  result = procrunner.run_process([
      'dials.export',
      'format=json',
      os.path.join(dials_regression, "centroid_test_data", "experiments.json"),
      os.path.join(dials_regression, "centroid_test_data", "integrated.pickle"),
      'json.filename=integrated.json',
      'n_digits=4',
      'compact=False',
  ])
  assert result['exitcode'] == 0
  assert result['stderr'] == ''
  assert os.path.exists('integrated.json')

  with open('integrated.json', 'rb') as f:
    d = json.load(f)
  assert "imageset_id" in d.keys()
  assert "rlp" in d.keys()
  assert "experiment_id" in d.keys()
  assert d['rlp'][:3] == [-0.5975, -0.6141, 0.4702], d['rlp'][:3]
  assert d['imageset_id'][0] == 0
  assert d['experiment_id'][0] == 0

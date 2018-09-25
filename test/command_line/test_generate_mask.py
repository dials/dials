from __future__ import absolute_import, division, print_function

import six.moves.cPickle as pickle
import os
import procrunner
import pytest

@pytest.fixture
def input_filename(dials_regression, run_in_tmpdir):
  yield os.path.join(dials_regression, "centroid_test_data", "experiments.json")
  print("temporary directory=", run_in_tmpdir.strpath)

def test_generate_mask(input_filename):
  result = procrunner.run_process([
      'dials.generate_mask',
      input_filename,
  ])
  assert result['exitcode'] == 0
  assert result['stderr'] == ''
  assert os.path.exists("mask.pickle")

def test_generate_mask_with_untrusted_rectangle(input_filename):
  result = procrunner.run_process([
      'dials.generate_mask',
      input_filename,
      'output.mask=mask2.pickle',
      'output.experiments=masked_experiments.json',
      'untrusted.rectangle=100,200,100,200'
  ])
  assert result['exitcode'] == 0
  assert result['stderr'] == ''
  assert os.path.exists("mask2.pickle")
  assert os.path.exists("masked_experiments.json")
  from dxtbx.serialize import load
  experiments = load.experiment_list("masked_experiments.json")
  imageset = experiments.imagesets()[0]
  assert imageset.external_lookup.mask.filename == os.path.join(
      os.path.abspath(os.getcwd()), 'mask2.pickle')

def test_generate_mask_with_untrusted_circle(input_filename):
  result = procrunner.run_process([
      'dials.generate_mask',
      input_filename,
      'output.mask=mask3.pickle',
      'untrusted.circle=100,100,10'
  ])
  assert result['exitcode'] == 0
  assert result['stderr'] == ''
  assert os.path.exists("mask3.pickle")

def test_generate_mask_with_resolution_range(input_filename):
  result = procrunner.run_process([
      'dials.generate_mask',
      input_filename,
      'output.mask=mask4.pickle',
      'resolution_range=2,3',
  ])
  assert result['exitcode'] == 0
  assert result['stderr'] == ''
  assert os.path.exists("mask4.pickle")

def test_generate_mask_with_d_min_d_max(input_filename):
  result = procrunner.run_process([
      'dials.generate_mask',
      input_filename,
      'output.mask=mask5.pickle',
      'd_min=3',
      'd_max=2',
  ])
  assert result['exitcode'] == 0
  assert result['stderr'] == ''
  assert os.path.exists("mask5.pickle")

def test_generate_mask_with_ice_rings(input_filename):
  result = procrunner.run_process([
      'dials.generate_mask',
      input_filename,
      'output.mask=mask6.pickle',
      'ice_rings{filter=True;d_min=2}',
  ])
  assert result['exitcode'] == 0
  assert result['stderr'] == ''
  assert os.path.exists("mask6.pickle")

def test_generate_mask_with_untrusted_polygon_and_pixels(input_filename):
  result = procrunner.run_process([
      'dials.generate_mask',
      input_filename,
      'output.mask=mask3.pickle',
      'untrusted.polygon=100,100,100,200,200,200,200,100',
      'untrusted.pixel=0,0',
      'untrusted.pixel=1,1'
  ])
  assert result['exitcode'] == 0
  assert result['stderr'] == ''
  assert os.path.exists("mask3.pickle")
  with open("mask3.pickle", "rb") as fh:
    mask = pickle.load(fh)
  assert not mask[0][0,0]
  assert not mask[0][1,1]
  assert mask[0][0,1]

from __future__ import absolute_import, division, print_function

from glob import glob
import os
from procrunner import run_process

def test_combining_spots(dials_regression, run_in_tmpdir):
  images = sorted(glob(os.path.join(dials_regression, 'centroid_test_data', "centroid*.cbf")))
  images_1 = images[0:int(len(images)/2)]
  images_2 = images[int(len(images)/2):]

  result = run_process(['dials.import'] + images_1 + ['output.experiments=experiments-1.json'])
  assert result['exitcode'] == 0
  assert result['stderr'] == ''
  assert os.path.exists('experiments-1.json')

  result = run_process(['dials.find_spots', 'experiments-1.json', 'output.reflections=strong-1.pickle'])
  assert result['exitcode'] == 0
  assert result['stderr'] == ''
  assert os.path.exists('strong-1.pickle')

  result = run_process(['dials.import'] + images_2 + ['output.experiments=experiments-2.json'])
  assert result['exitcode'] == 0
  assert result['stderr'] == ''
  assert os.path.exists('experiments-2.json')

  result = run_process(['dials.find_spots', 'experiments-2.json', 'output.reflections=strong-2.pickle'])
  assert result['exitcode'] == 0
  assert result['stderr'] == ''
  assert os.path.exists('strong-2.pickle')

  result = run_process([
      'dials.combine_found_spots',
      'experiments-1.json',
      'experiments-2.json',
      'strong-1.pickle',
      'strong-2.pickle',
      'output.reflections=combined.pickle',
      'output.experiments=combined.json',
  ])
  assert result['exitcode'] == 0
  assert result['stderr'] == ''
  assert os.path.exists('combined.json')
  assert os.path.exists('combined.pickle')

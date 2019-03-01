from __future__ import absolute_import, division, print_function

from glob import glob
import os
from procrunner import run_process
from dials.array_family import flex

def test_combining_spots(dials_regression, run_in_tmpdir):
  images = sorted(glob(os.path.join(dials_regression, 'centroid_test_data', "centroid*.cbf")))
  images_1 = images[0:int(len(images)/2)]
  images_2 = images[int(len(images)/2):]

  result = run_process(['dials.import'] + images_1 + ['output.datablock=datablock-1.json'])
  assert result['exitcode'] == 0
  assert result['stderr'] == ''
  assert os.path.exists('datablock-1.json')

  result = run_process(['dials.find_spots', 'datablock-1.json', 'output.reflections=strong-1.pickle'])
  assert result['exitcode'] == 0
  assert result['stderr'] == ''
  assert os.path.exists('strong-1.pickle')

  result = run_process(['dials.import'] + images_2 + ['output.datablock=datablock-2.json'])
  assert result['exitcode'] == 0
  assert result['stderr'] == ''
  assert os.path.exists('datablock-2.json')

  result = run_process(['dials.find_spots', 'datablock-2.json', 'output.reflections=strong-2.pickle'])
  assert result['exitcode'] == 0
  assert result['stderr'] == ''
  assert os.path.exists('strong-2.pickle')

  result = run_process([
      'dials.combine_found_spots',
      'datablock-1.json',
      'datablock-2.json',
      'strong-1.pickle',
      'strong-2.pickle',
      'output.reflections=combined.pickle',
      'output.datablock=combined.json',
  ])
  assert result['exitcode'] == 0
  assert result['stderr'] == ''
  assert os.path.exists('combined.json')
  assert os.path.exists('combined.pickle')

  r = flex.reflection_table.from_pickle('combined.pickle')
  assert r['id'].all_eq(0)

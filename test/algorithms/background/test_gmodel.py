from __future__ import absolute_import, division, print_function

def test_run(dials_regression, run_in_tmpdir):
  import os
  path = os.path.join(dials_regression, "centroid_test_data")

  from libtbx import easy_run

  os.mkdir('simple')
  os.mkdir('robust')
  os.mkdir('gmodel_simple')
  os.mkdir('gmodel_robust')

  assert os.path.exists(os.path.join(path, 'experiments.json'))

  from dials.array_family import flex
  from dials.algorithms.background.gmodel import StaticBackgroundModel
  ysize = 2527
  xsize = 2463
  data = flex.double(flex.grid(ysize, xsize), 1)
  model = StaticBackgroundModel()
  model.add(data)
  import six.moves.cPickle as pickle
  with open("model.pickle", "wb") as fh:
    pickle.dump(model, fh, pickle.HIGHEST_PROTOCOL)

  # Call dials.integrate
  easy_run.fully_buffered([
    'dials.integrate',
    os.path.join(path, 'experiments.json'),
    'profile.fitting=False',
    'background.algorithm=simple',
    'background.simple.outlier.algorithm=null',
    'output.reflections=./simple/reflections.pickle'
  ]).raise_if_errors()

  assert os.path.exists("simple/reflections.pickle")

  # Call dials.integrate
  easy_run.fully_buffered([
    'dials.integrate',
    os.path.join(path, 'experiments.json'),
    'profile.fitting=False',
    'background.algorithm=glm',
    'output.reflections=./robust/reflections.pickle'
  ]).raise_if_errors()

  assert os.path.exists("robust/reflections.pickle")

  # Call dials.integrate
  easy_run.fully_buffered([
    'dials.integrate',
    os.path.join(path, 'experiments.json'),
    'profile.fitting=False',
    'background.algorithm=gmodel',
    'background.gmodel.robust.algorithm=False',
    'background.gmodel.model=model.pickle',
    'output.reflections=./gmodel_simple/reflections.pickle'
  ]).raise_if_errors()

  assert os.path.exists("gmodel_simple/reflections.pickle")

  # Call dials.integrate
  easy_run.fully_buffered([
    'dials.integrate',
    os.path.join(path, 'experiments.json'),
    'profile.fitting=False',
    'background.algorithm=gmodel',
    'background.gmodel.robust.algorithm=True',
    'background.gmodel.model=model.pickle',
    'output.reflections=./gmodel_robust/reflections.pickle'
  ]).raise_if_errors()

  assert os.path.exists("gmodel_robust/reflections.pickle")

  reflections1 = flex.reflection_table.from_pickle("simple/reflections.pickle")
  reflections2 = flex.reflection_table.from_pickle("robust/reflections.pickle")
  reflections3 = flex.reflection_table.from_pickle("gmodel_simple/reflections.pickle")
  reflections4 = flex.reflection_table.from_pickle("gmodel_robust/reflections.pickle")

  assert len(reflections1) == len(reflections3)
  assert len(reflections2) == len(reflections4)

  flag = flex.reflection_table.flags.integrated_sum
  integrated1 = reflections1.select(reflections1.get_flags(flag, all=True))
  integrated2 = reflections2.select(reflections2.get_flags(flag, all=True))
  integrated3 = reflections3.select(reflections3.get_flags(flag, all=True))
  integrated4 = reflections4.select(reflections4.get_flags(flag, all=True))

  assert len(integrated1) > 0
  assert len(integrated2) > 0
  assert len(integrated1) == len(integrated3)
  assert len(integrated2) == len(integrated4)

  mean_bg1 = integrated1['background.mean']
  mean_bg2 = integrated2['background.mean']
  mean_bg3 = integrated3['background.mean']
  mean_bg4 = integrated4['background.mean']

  scale3 = integrated3['background.scale']
  scale4 = integrated4['background.scale']

  diff1 = flex.abs(mean_bg1 - mean_bg3)
  diff2 = flex.abs(mean_bg2 - mean_bg4)
  assert (scale3 > 0).count(False) == 0
  assert (scale4 > 0).count(False) == 0
  assert (diff1 < 1e-5).count(False) == 0

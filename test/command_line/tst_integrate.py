from __future__ import division
from dials.array_family import flex # import dependency


class Test(object):

  def __init__(self):
    from os.path import join
    import libtbx.load_env
    try:
      dials_regression = libtbx.env.dist_path('dials_regression')
    except KeyError, e:
      print 'SKIP: dials_regression not configured'
      exit(0)

    self.path = join(dials_regression, "centroid_test_data")

  def run(self):
    # self.test1()
    # TODO: Why is this test disabled? Fix or remove?
    self.test2()
    self.test3()
    self.test4()

  def test1(self):
    from os.path import join, exists
    from libtbx import easy_run
    import os
    from uuid import uuid4

    dirname ='tmp_%s' % uuid4().hex
    os.mkdir(dirname)
    os.chdir(dirname)

    assert exists(join(self.path, 'experiments.json'))

    # Call dials.integrate
    easy_run.fully_buffered([
      'dials.integrate',
      join(self.path, 'experiments.json'),
      'profile.fitting=False',
    ]).raise_if_errors()

    import cPickle as pickle
    table = pickle.load(open('integrated.pickle', 'rb'))
    print len(table)
    assert(len(table) == 751)

    # assert(len(table) == 764)
    assert('id' in table)
    for row in table:
      assert(row['id'] == 0)
    self.table = table
    print 'OK'

  def test2(self):
    from os.path import join
    from libtbx import easy_run
    import os

    # Call dials.integrate
    easy_run.fully_buffered([
      'dials.integrate',
      join(self.path, 'experiments.json'),
      'profile.fitting=False',
      'integration.integrator=3d',
    ]).raise_if_errors()

    import cPickle as pickle
    table = pickle.load(open('integrated.pickle', 'rb'))
    mask = table.get_flags(table.flags.integrated,all=False)
    assert(len(table) == 1996)
    assert(mask.count(True) == 1666)

    # assert(len(table) == 764)
    assert('id' in table)
    for row in table:
      assert(row['id'] == 0)
    self.table = table
    print 'OK'

  def test3(self):
    from os.path import join
    from libtbx import easy_run
    import shutil
    from os.path import join
    shutil.copyfile(join(self.path, "experiments.json"), "experiments.json")
    for i in range(1, 10):
      shutil.copyfile(join(self.path, "centroid_000%d.cbf" % i),
                      "centroid_001%d.cbf" % i)

    infile = open("experiments.json", "r")
    lines = infile.readlines()
    inscan = False
    inimagerange = False
    inoscillation = False
    count = 0
    done1 = False
    done2 = False
    for i, line in enumerate(lines):
      if not inscan:
        if line.strip().startswith('"scan": ['):
          inscan = True
      else:
        if not inimagerange and not inoscillation:
          if line.strip().startswith('"image_range": ['):
            inimagerange = True
          if line.strip().startswith('"oscillation": ['):
            inoscillation = True
        elif inimagerange:
          if count == 0:
            lines[i] = '11,'
            count += 1
          elif count == 1:
            lines[i] = '19'
            done1 = True
            inimagerange = False
            count = 0
        elif inoscillation:
          if count == 0:
            lines[i] = '360.0,'
            done2 = True
            inoscillation = False
            inscan = False
            break
    assert(done1 == True)
    assert(done2 == True)
    infile.close()
    outfile = open("experiments.json", "w")
    outfile.write('\n'.join(lines))
    outfile.close()

    # Call dials.integrate
    easy_run.fully_buffered([
      'dials.integrate',
      'experiments.json',
      'profile.fitting=False',
      'integration.integrator=3d',
    ]).raise_if_errors()

    from math import pi
    import cPickle as pickle
    table = pickle.load(open('integrated.pickle', 'rb'))
    mask1 = table.get_flags(table.flags.integrated,all=False)
    assert(len(table) == 1996), "Table has %d entries instead of 1996" % len(table)
    assert(mask1.count(True) == 1666)
    mask2 = self.table.get_flags(table.flags.integrated,all=False)
    assert(mask1.all_eq(mask2))
    t1 = table.select(mask1)
    t2 = self.table.select(mask1)
    Cal_P1 = t1['xyzcal.mm'].parts()[2]
    Cal_Z1 = t1['xyzcal.px'].parts()[2]
    Obs_Z1 = t1['xyzobs.px.value'].parts()[2]
    # Obs_P1 = t1['xyzobs.mm.value'].parts()[2]
    Cal_Z2 = t2['xyzcal.px'].parts()[2]
    Cal_P2 = t2['xyzcal.mm'].parts()[2]
    Obs_Z2 = t2['xyzobs.px.value'].parts()[2]
    # Obs_P2 = t2['xyzobs.mm.value'].parts()[2]
    diff_I = t1['intensity.sum.value'] - t2['intensity.sum.value']
    diff_Cal_Z = Cal_Z1 - (Cal_Z2 + 10)
    diff_Obs_Z = Obs_Z1 - (Obs_Z2 + 10)
    diff_Cal_P = Cal_P1 - (Cal_P2 + 2*pi)
    # diff_Obs_P = Obs_P1 - (Obs_P2 + 2*pi)
    assert(flex.abs(diff_I).all_lt(1e-7))
    assert(flex.abs(diff_Cal_Z).all_lt(1e-7))
    assert(flex.abs(diff_Cal_P).all_lt(1e-7))
    assert(flex.abs(diff_Obs_Z).all_lt(1e-7))
    # assert(flex.abs(diff_Obs_P).all_lt(1e-7))

    print 'OK'

  def test4(self):
    from os.path import join
    from libtbx import easy_run
    import os

    dirname ='test4'
    os.mkdir(dirname)
    os.chdir(dirname)

    # Call dials.integrate
    easy_run.fully_buffered([
      'dials.integrate',
      join(self.path, 'experiments.json'),
      'profile.fitting=False',
      'sampling.integrate_all_reflections=False',
    ]).raise_if_errors()

    import cPickle as pickle
    table = pickle.load(open('integrated.pickle', 'rb'))
    assert len(table) == 1000

    # Call dials.integrate
    easy_run.fully_buffered([
      'dials.integrate',
      join(self.path, 'experiments.json'),
      'profile.fitting=False',
      'sampling.integrate_all_reflections=False',
      'sampling.minimum_sample_size=500',
    ]).raise_if_errors()

    import cPickle as pickle
    table = pickle.load(open('integrated.pickle', 'rb'))
    assert len(table) == 500

    print 'OK'

if __name__ == '__main__':
  from dials.test import cd_auto
  with cd_auto(__file__):
    test = Test()
    test.run()

from __future__ import absolute_import, division
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
    self.integration_test_data = join(dials_regression, "integration_test_data")

  def run(self):
    # self.test1()
    # TODO: Why is this test disabled? Fix or remove?
    self.test2()
    self.test3()
    self.test4()
    self.test_multi_sweep()
    self.test_multi_lattice()
    self.test_output_rubbish()

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

    with open("experiments.json", "r") as infile:
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
    assert(done1)
    assert(done2)
    with open("experiments.json", "w") as outfile:
      outfile.write('\n'.join(lines))

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

  def test_multi_sweep(self):
    from os.path import join
    from libtbx import easy_run
    import os

    dirname ='multi_sweep'
    os.mkdir(dirname)
    os.chdir(dirname)

    # Call dials.integrate
    easy_run.fully_buffered([
      'dials.integrate',
      join(self.integration_test_data, 'multi_sweep', 'experiments.json'),
      join(self.integration_test_data, 'multi_sweep', 'indexed.pickle')
    ]).raise_if_errors()

    import cPickle as pickle
    table = pickle.load(open('integrated.pickle', 'rb'))
    assert len(table) == 4020

    # Check the results
    T1 = table[:2010]
    T2 = table[2010:]
    ID1 = list(set(T1['id']))
    ID2 = list(set(T2['id']))
    assert len(ID1) == 1
    assert len(ID2) == 1
    assert ID1[0] == 0
    assert ID2[0] == 1
    I1 = T1['intensity.prf.value']
    I2 = T2['intensity.prf.value']
    F1 = T1.get_flags(T1.flags.integrated_prf)
    F2 = T2.get_flags(T2.flags.integrated_prf)
    assert F1 == F2
    I1 = I1.select(F1)
    I2 = I2.select(F2)
    assert flex.abs(I1 - I2) < 1e-6

    print 'OK'

  def test_multi_lattice(self):
    from os.path import join
    from libtbx import easy_run
    import os

    dirname ='multi_sweep'
    os.mkdir(dirname)
    os.chdir(dirname)

    # Call dials.integrate
    easy_run.fully_buffered([
      'dials.integrate',
      join(self.integration_test_data, 'multi_lattice', 'experiments.json'),
      join(self.integration_test_data, 'multi_lattice', 'indexed.pickle')
    ]).raise_if_errors()

    import cPickle as pickle
    table = pickle.load(open('integrated.pickle', 'rb'))
    assert len(table) == 5605, "%d, %d" % (len(table), 5605)

    # Check output contains from two lattices
    exp_id = list(set(table['id']))
    assert len(exp_id) == 2

    # Check both lattices have integrated reflections
    mask = table.get_flags(table.flags.integrated_prf)
    table = table.select(mask)
    exp_id = list(set(table['id']))
    assert len(exp_id) == 2

    print 'OK'


  def test_output_rubbish(self):
    from os.path import join, exists
    from libtbx import easy_run
    import os
    from uuid import uuid4

    dirname ='tmp_%s' % uuid4().hex
    os.mkdir(dirname)
    os.chdir(dirname)

    assert exists(join(self.path, 'datablock.json'))
    assert exists(join(self.path, 'strong.pickle'))

    # Call dials.integrate
    easy_run.fully_buffered([
      'dials.index',
      join(self.path, 'datablock.json'),
      join(self.path, 'strong.pickle'),
    ]).raise_if_errors()

    assert exists('experiments.json')
    assert exists('indexed.pickle')

    # Call dials.integrate
    easy_run.fully_buffered([
      'dials.integrate',
      'experiments.json',
      'indexed.pickle',
      'profile.fitting=False',
    ]).raise_if_errors()

    import cPickle as pickle
    table = pickle.load(open('integrated.pickle', 'rb'))
    assert table.get_flags(table.flags.bad_reference) > 0

    assert('id' in table)
    for row in table:
      assert(row['id'] == 0)
    self.table = table
    print 'OK'


if __name__ == '__main__':
  from dials.test import cd_auto
  with cd_auto(__file__):
    test = Test()
    test.run()

from __future__ import division


class Test(object):

  def __init__(self):
    from os.path import join
    from dials.array_family import flex
    import libtbx.load_env
    try:
      dials_regression = libtbx.env.dist_path('dials_regression')
    except KeyError, e:
      print 'FAIL: dials_regression not configured'
      exit(0)

    self.path = join(dials_regression, "centroid_test_data")

  def run(self):
    self.test1()
    self.test2()
    self.test3()

  def test1(self):
    from os.path import abspath, join, exists
    from libtbx import easy_run
    import os
    from uuid import uuid4

    dirname ='tmp_%s' % uuid4().hex
    os.mkdir(dirname)
    os.chdir(dirname)

    assert exists(join(self.path, 'experiments.json'))
    assert exists(join(self.path, 'profile.phil'))

    # Call dials.integrate
    easy_run.fully_buffered([
      'dials.integrate',
      join(self.path, 'experiments.json'),
      join(self.path, 'profile.phil'),
      'intensity.algorithm=sum3d',
    ]).raise_if_errors()

    import cPickle as pickle
    table = pickle.load(open('integrated.pickle', 'rb'))
    mask = table.get_flags(table.flags.integrated,all=False)
    assert(len(table) == 1996)
    assert(mask.count(True) == 1995)

    # assert(len(table) == 764)
    assert('id' in table)
    for row in table:
      assert(row['id'] == 0)
    self.table = table
    print 'OK'

  def test2(self):
    from os.path import abspath, join
    from libtbx import easy_run
    import os
    from uuid import uuid4

    dirname ='tmp_%s' % uuid4().hex
    os.mkdir(dirname)
    os.chdir(dirname)

    # Call dials.integrate
    easy_run.fully_buffered([
      'dials.integrate2',
      join(self.path, 'experiments.json'),
      join(self.path, 'profile.phil'),
      'intensity.algorithm=sum3d',
    ]).raise_if_errors()

    import cPickle as pickle
    table = pickle.load(open('integrated.pickle', 'rb'))
    mask = table.get_flags(table.flags.integrated,all=False)
    assert(len(table) == 1996)
    assert(mask.count(True) == 1995)

    # assert(len(table) == 764)
    assert('id' in table)
    for row in table:
      assert(row['id'] == 0)
    self.table = table
    print 'OK'

  def test3(self):
    from os.path import abspath, join
    from libtbx import easy_run
    import os
    from uuid import uuid4
    import shutil
    from os.path import join
    from dials.array_family import flex
    shutil.copyfile(join(self.path, "experiments.json"), "experiments.json")
    shutil.copyfile(join(self.path, "profile.phil"), "profile.phil")
    for i in range(1, 10):
      shutil.copyfile(join(self.path, "centroid_000%d.cbf" % i),
                      "centroid_001%d.cbf" % i)

    infile = open("experiments.json", "r")
    lines = infile.readlines()
    inscan = False
    inimagerange = False
    count = 0
    done = False
    for i, line in enumerate(lines):
      if not inscan:
        if line.strip().startswith('"scan": ['):
          inscan = True
      else:
        if not inimagerange:
          if line.strip().startswith('"image_range": ['):
            inimagerange = True
        else:
          if count == 0:
            lines[i] = '11,'
            count += 1
          elif count == 1:
            lines[i] = '19'
            done = True
            break
    assert(done == True)
    infile.close()
    outfile = open("experiments.json", "w")
    outfile.write('\n'.join(lines))
    outfile.close()

    # Call dials.integrate
    easy_run.fully_buffered([
      'dials.integrate2',
      'experiments.json',
      'profile.phil',
      'intensity.algorithm=sum3d',
    ]).raise_if_errors()

    import cPickle as pickle
    table = pickle.load(open('integrated.pickle', 'rb'))
    mask = table.get_flags(table.flags.integrated,all=False)
    assert(len(table) == 1996)
    assert(mask.count(True) == 1995)
    Cal_Z1 = table['xyzcal.px'].parts()[2]
    Cal_Z2 = self.table['xyzcal.px'].parts()[2]
    Obs_Z1 = table['xyzobs.px.value'].parts()[2]
    Obs_Z2 = self.table['xyzobs.px.value'].parts()[2]
    diff_I = table['intensity.sum.value'] - self.table['intensity.sum.value']
    diff_Cal_Z = Cal_Z1 - Cal_Z2
    diff_Obs_Z = Obs_Z1 - Obs_Z2
    assert((flex.abs(diff_I)     - 10).all_lt(1e-7))
    assert((flex.abs(diff_Cal_Z) - 10).all_lt(1e-7))
    assert((flex.abs(diff_Obs_Z) - 10).all_lt(1e-7))

    print 'OK'

if __name__ == '__main__':
  from dials.test import cd_auto
  with cd_auto(__file__):
    test = Test()
    test.run()

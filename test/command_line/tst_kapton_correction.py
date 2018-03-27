from __future__ import absolute_import, division

from dials.array_family import flex  # import dependency

class Test(object):

  def __init__(self):
    from os.path import join
    import libtbx.load_env
    try:
      dials_regression = libtbx.env.dist_path('dials_regression')
    except KeyError:
      print 'SKIP: dials_regression not configured'
      exit(0)

    self.path = join(dials_regression, "integration_test_data", "stills_PSII")

  def run(self):
    self.test_integrate_with_kapton()

  def test_integrate_with_kapton(self):
    from os.path import join, exists
    from libtbx import easy_run
    import os
    from uuid import uuid4

    dirname ='tmp_%s' % uuid4().hex
    os.mkdir(dirname)
    os.chdir(dirname)
    loc = os.getcwd()

    pickle_name = 'idx-20161021225550223_indexed.pickle'
    json_name = 'idx-20161021225550223_refined_experiments.json'
    image_name = '20161021225550223.pickle'
    pickle_path = join(self.path, pickle_name)
    json_path = join(self.path, json_name)
    image_path = join(self.path, image_name)
    assert os.path.exists(pickle_path)
    assert os.path.exists(json_path)

    import shutil
    shutil.copy(pickle_path, loc)
    shutil.copy(image_path, loc)

    w = open(json_name, 'wb')
    r = open(json_path, 'rb')
    w.write(r.read() % loc)
    r.close()
    w.close()

    templ_phil = """
      output {
        experiments = 'idx-20161021225550223_integrated_experiments_%s.json'
        reflections = 'idx-20161021225550223_integrated_%s.pickle'
      }
      integration {
        lookup.mask = '%s'
        integrator = stills
        profile.fitting = False
        background.algorithm = simple
        debug {
          output = True
          separate_files = False
          split_experiments = False
        }
      }
      profile {
        gaussian_rs.min_spots.overall = 0
      }
      absorption_correction {
        apply = %s
        algorithm = fuller_kapton
        fuller_kapton {
          smart_sigmas = True
        }
      }
"""
    without_kapton_phil = templ_phil % ("nokapton", "nokapton",
      os.path.join(self.path, "mask.pickle"), "False")
    with_kapton_phil = templ_phil % ("kapton", "kapton",
      os.path.join(self.path, "mask.pickle"), "True")

    f = open("integrate_without_kapton.phil", 'wb')
    f.write(without_kapton_phil)
    f.close()

    f = open("integrate_with_kapton.phil", 'wb')
    f.write(with_kapton_phil)
    f.close()



    # Call dials.integrate with and without kapton correction
    for phil in "integrate_without_kapton.phil", "integrate_with_kapton.phil":
      result = easy_run.fully_buffered([
        'dials.integrate', pickle_name, json_name, phil
      ]).raise_if_errors()

    import cPickle as pickle
    results = []
    for mode in "kapton", "nokapton":
      result = os.path.join(loc, "idx-20161021225550223_integrated_%s.pickle" % mode)
      table = pickle.load(open(result, 'rb'))
      millers = table['miller_index']
      test_indices = {'zero':(-5, 2, -6), 'low':(-2, -20, 7), 'high':(-1, -10, 4)}
      test_rows = {k:millers.first_index(v) for k,v in test_indices.iteritems()}
      test_I_sigsqI = {k:(table[v]['intensity.sum.value'], table[v]['intensity.sum.variance'])
                       for k,v in test_rows.iteritems()}
      results.append(test_I_sigsqI)
    assert results[0]['zero'][0] == results[1]['zero'][0]
    assert results[0]['zero'][1] - results[1]['zero'][1] < 0.0001
    assert False not in [results[0]['low'][i] > results[1]['low'][i] for i in (0, 1)]
    assert False not in [results[0]['high'][i] > results[1]['high'][i] for i in (0, 1)]

if __name__ == '__main__':
  from dials.test import cd_auto
  with cd_auto(__file__):
    test = Test()
    test.run()

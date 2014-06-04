from __future__ import division


class TestImporter:

  def __init__(self):

    import os
    import libtbx.load_env
    try:
      dials_regression = libtbx.env.dist_path('dials_regression')
    except KeyError, e:
      print 'FAIL: dials_regression not configured'
      exit(0)

    self.path = os.path.join(dials_regression, 'centroid_test_data')
    self.create_data()

  def create_data(self):
    from dials.array_family import flex
    from dials.model.data import ReflectionList
    import cPickle as pickle
    table = flex.reflection_table()
    table['col1'] = flex.int(10)
    table['col2'] = flex.int(10)
    pickle.dump(table, open('test_reflections1.p', 'wb'))

    table = flex.reflection_table()
    table['col2'] = flex.int(10)
    table['col3'] = flex.int(10)
    pickle.dump(table, open('test_reflections2.p', 'wb'))


  def run(self):
    from glob import glob
    import os

    arguments = glob(os.path.join(self.path, 'centroid*.cbf'))

    arguments = list(arguments) + [
        os.path.join(self.path, 'datablock.json'),
        os.path.join(self.path, 'experiments.json'),
        os.path.join(self.path, 'non_existent.file'),
        os.path.join(self.path, 'datablock.json'),
        'test_reflections1.p',
        'test_reflections2.p',
        os.path.join(self.path, 'another_non_existent.file'),
        os.path.join(self.path, 'experiments.json'),
    ]

    from dials.util.command_line import Importer
    importer = Importer(arguments, verbose=False)
    assert(importer.reflections is not None)
    assert(len(importer.reflections) == 2)
    assert(importer.reflections[0].ncols() == 2)
    assert(importer.reflections[1].ncols() == 2)
    assert(importer.reflections[0].nrows() == 10)
    assert(importer.reflections[1].nrows() == 10)
    assert(len(importer.experiments) == 2)
    assert(len(importer.datablocks) == 3)
    assert(len(importer.unhandled_arguments) == 2)

    print 'OK'

if __name__ == '__main__':
  from dials.test import cd_auto
  with cd_auto(__file__):
    test = TestImporter()
    test.run()

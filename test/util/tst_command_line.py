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

    arguments = [
        os.path.join(self.path, 'sweep.json'),
        os.path.join(self.path, 'crystal.json'),
        os.path.join(self.path, 'non_existent.file'),
        os.path.join(self.path, 'crystal.json'),
        'test_reflections1.p',
        'test_reflections2.p',
        os.path.join(self.path, 'extracted.tar'),
        os.path.join(self.path, 'another_non_existent.file'),
        os.path.join(self.path, 'sweep.json'),
    ]
    arguments.extend(glob(os.path.join(self.path, 'centroid_*.cbf')))

    from dials.util.command_line import Importer
    importer = Importer(arguments, verbose=False)

    assert(importer.reflections.ncols() == 3)
    assert(importer.reflections.nrows() == 10)
    assert(len(importer.imagesets) == 3)
    assert(len(importer.crystals) == 2)

    if os.path.exists(os.path.join(self.path, 'extracted.tar')):
      assert(importer.extracted != None)
      assert(len(importer.unhandled) == 2)
    else:
      assert(len(importer.unhandled) == 3)

    print 'OK'

if __name__ == '__main__':
  test = TestImporter()
  test.run()

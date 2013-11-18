

class TestImporter:

  def __init__(self):

    import os
    import libtbx.load_env
    try:
      dials_regression = libtbx.env.dist_path('dials_regression')
    except KeyError, e:
      print 'FAIL: dials_regression not configured'
      return

    self.path = os.path.join(dials_regression, 'centroid_test_data')

  def run(self):
    from glob import glob
    import os

    arguments = [
        os.path.join(self.path, 'sweep.json'),
        os.path.join(self.path, 'crystal.json'),
        os.path.join(self.path, 'non_existent.file'),
        os.path.join(self.path, 'crystal.json'),
        os.path.join(self.path, 'spot_xds.pickle'),
        os.path.join(self.path, 'extracted.tar'),
        os.path.join(self.path, 'another_non_existent.file'),
        os.path.join(self.path, 'sweep.json'),
    ]
    arguments.extend(glob(os.path.join(self.path, 'centroid_*.cbf')))

    from dials.util.command_line import Importer
    importer = Importer(arguments)
    assert(len(importer.reflections) == 664)
    assert(len(importer.imagesets) == 3)
    assert(len(importer.crystals) == 2)

    if os.path.exists(os.path.join(self.path, 'extracted.tar')):
      assert(importer.extracted != None)
      assert(len(importer.unhandled_arguments) == 2)
    else:
      assert(len(importer.unhandled_arguments) == 3)

    print 'OK'

if __name__ == '__main__':
  test = TestImporter()
  test.run()

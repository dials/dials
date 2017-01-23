
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

  def run(self):

    from os.path import join
    from libtbx import easy_run

    input_filename = join(self.path, "datablock.json")
    output_filename = "output_datablock.json"
    mask_filename = join(self.path, "lookup_mask.pickle")

    easy_run.fully_buffered(
      ['dials.apply_mask',
       'input.datablock=%s' % input_filename,
       'input.mask=%s' % mask_filename,
       'output.datablock=%s' % output_filename]).raise_if_errors()

    from dxtbx.datablock import DataBlockFactory
    datablocks = DataBlockFactory.from_json_file(output_filename)

    assert len(datablocks) == 1
    imagesets = datablocks[0].extract_imagesets()
    assert len(imagesets) == 1
    imageset = imagesets[0]
    assert imageset.external_lookup.mask.filename == mask_filename

    print 'OK'

if __name__ == '__main__':
  from dials.test import cd_auto
  with cd_auto(__file__):
    test = Test()
    test.run()

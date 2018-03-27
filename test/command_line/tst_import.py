from __future__ import absolute_import, division

from libtbx.test_utils import approx_equal

class Test(object):
  def __init__(self):
    from os.path import join
    import libtbx.load_env
    try:
      dials_regression = libtbx.env.dist_path('dials_regression')
    except KeyError:
      print 'FAIL: dials_regression not configured'
      exit(0)

    self.path = join(dials_regression, "centroid_test_data")

  def run(self):
    self.tst_from_image_files()
    self.tst_import_beam_centre()
    self.tst_with_mask()
    self.tst_override_geometry()
    self.tst_extrapolate_scan()
    self.tst_multiple_sweeps()

  def tst_multiple_sweeps(self):
    from glob import glob
    import os
    from libtbx import easy_run
    from dxtbx.serialize import load

    # Find the image files
    image_files = sorted(glob(os.path.join(self.path, "centroid*.cbf")))
    del image_files[4] # Delete filename to force two sweeps
    image_files = ' '.join(image_files)

    # provide mosflm beam centre to dials.import
    cmd = 'dials.import %s output.datablock=datablock_multiple_sweeps.json' % (image_files)
    try:
      easy_run.fully_buffered(cmd).raise_if_errors()
      assert False, "Expected exception"
    except Exception:
      pass

    cmd = 'dials.import %s output.datablock=datablock_multiple_sweeps.json allow_multiple_sweeps=True' % (image_files)
    easy_run.fully_buffered(cmd).raise_if_errors()
    assert os.path.exists("datablock_multiple_sweeps.json")
    datablock = load.datablock("datablock_multiple_sweeps.json")[0]
    imgset = datablock.extract_imagesets()
    assert len(imgset) == 2


  def tst_with_mask(self):
    from glob import glob
    import os
    from libtbx import easy_run
    from dxtbx.serialize import load

    # Find the image files
    image_files = glob(os.path.join(self.path, "centroid*.cbf"))
    image_files = ' '.join(image_files)
    mask_filename = os.path.join(self.path, "mask.pickle")

    # provide mosflm beam centre to dials.import
    cmd = 'dials.import %s mask=%s output.datablock=datablock_with_mask.json' % (image_files, mask_filename)
    easy_run.fully_buffered(cmd).raise_if_errors()
    assert os.path.exists("datablock_with_mask.json")
    datablock = load.datablock("datablock_with_mask.json")[0]
    imgset = datablock.extract_imagesets()[0]
    assert imgset.external_lookup.mask.filename == mask_filename


  def tst_override_geometry(self):
    from glob import glob
    import os
    from libtbx import easy_run
    from dxtbx.serialize import load

    # Find the image files
    image_files = glob(os.path.join(self.path, "centroid*.cbf"))
    image_files = ' '.join(image_files)

    # Write a geometry phil file
    with open("geometry.phil", "w") as outfile:
      outfile.write(
        '''
        geometry {
          beam {
            wavelength = 2
            direction = (-1,0,0)
          }
          detector {
            panel {
              name = "New panel"
              type = "New type"
              pixel_size = 10,20
              image_size = 30,40
              trusted_range = 50,60
              thickness = 70
              material = "Si"
              fast_axis = -1,0,0
              slow_axis = 0,-1,0
              origin = 100,100,100
            }
          }
          goniometer {
            axes = 0,0,-1
            fixed_rotation = 0,1,2,3,4,5,6,7,8
            setting_rotation = 8,7,6,5,4,3,2,1,0
          }
          scan {
            image_range = 1,4
            oscillation = 1,2
          }
        }
    ''')

    # provide mosflm beam centre to dials.import
    cmd = 'dials.import %s geometry.phil output.datablock=override_geometry.json' %image_files
    easy_run.fully_buffered(cmd).raise_if_errors()
    assert os.path.exists("override_geometry.json")
    datablock = load.datablock("override_geometry.json")[0]
    imgset = datablock.extract_imagesets()[0]

    beam = imgset.get_beam()
    detector = imgset.get_detector()
    goniometer = imgset.get_goniometer()
    scan = imgset.get_scan()

    assert beam.get_wavelength() == 2
    assert beam.get_direction() == (-1,0,0)
    assert detector[0].get_name() == "New panel"
    assert detector[0].get_type() == "New type"
    assert detector[0].get_pixel_size() == (10,20)
    assert detector[0].get_image_size() == (30,40)
    assert detector[0].get_trusted_range() == (50,60)
    assert detector[0].get_thickness() == 70
    assert detector[0].get_material() == "Si"
    assert detector[0].get_fast_axis() == (-1,0,0)
    assert detector[0].get_slow_axis() == (0,-1,0)
    assert detector[0].get_origin() == (100,100,100)
    assert goniometer.get_rotation_axis_datum() == (0,0,-1)
    assert goniometer.get_fixed_rotation() == (0,1,2,3,4,5,6,7,8)
    assert goniometer.get_setting_rotation() == (8,7,6,5,4,3,2,1,0)
    assert scan.get_image_range() == (1,4)
    assert scan.get_oscillation() == (1,2)


  def tst_import_beam_centre(self):
    from glob import glob
    import os
    from libtbx import easy_run
    from dxtbx.serialize import load

    # Find the image files
    image_files = glob(os.path.join(self.path, "centroid*.cbf"))
    image_files = ' '.join(image_files)

    # provide mosflm beam centre to dials.import
    cmd = 'dials.import %s mosflm_beam_centre=100,200 output.datablock=mosflm_beam_centre.json' %image_files
    easy_run.fully_buffered(cmd).raise_if_errors()
    assert os.path.exists("mosflm_beam_centre.json")
    datablock = load.datablock("mosflm_beam_centre.json")[0]
    imgset = datablock.extract_imagesets()[0]
    beam_centre = imgset.get_detector()[0].get_beam_centre(imgset.get_beam().get_s0())
    assert approx_equal(beam_centre, (200,100))

    # provide an alternative datablock.json to get geometry from
    cmd = 'dials.import %s reference_geometry=mosflm_beam_centre.json output.datablock=mosflm_beam_centre2.json' %image_files
    easy_run.fully_buffered(cmd).raise_if_errors()
    assert os.path.exists("mosflm_beam_centre2.json")
    datablock = load.datablock("mosflm_beam_centre2.json")[0]
    imgset = datablock.extract_imagesets()[0]
    beam_centre = imgset.get_detector()[0].get_beam_centre(imgset.get_beam().get_s0())
    assert approx_equal(beam_centre, (200,100))

    # test slow_fast_beam_centre with a multi-panel CS-PAD image
    dials_regression = os.path.join(self.path, os.pardir)
    impath = os.path.join(dials_regression, "image_examples",
        "LCLS_cspad_nexus", "idx-20130301060858401.cbf")
    cmd = ('dials.import {0} slow_fast_beam_centre=134,42,18 '
           'output.datablock=slow_fast_beam_centre.json').format(impath)
    easy_run.fully_buffered(cmd).raise_if_errors()
    assert os.path.exists('slow_fast_beam_centre.json')
    datablock = load.datablock('slow_fast_beam_centre.json')[0]
    imgset = datablock.extract_imagesets()[0]
    # beam centre on 18th panel
    s0 = imgset.get_beam().get_s0()
    beam_centre = imgset.get_detector()[18].get_beam_centre_px(s0)
    assert approx_equal(beam_centre, (42,134))
    # check relative panel positions have not changed
    from scitbx import matrix
    o = matrix.col(imgset.get_detector()[0].get_origin())
    offsets = []
    for p in imgset.get_detector():
      intra_pnl = o - matrix.col(p.get_origin())
      offsets.append(intra_pnl.length())
    cmd = ('dials.import {0} output.datablock=reference.json').format(impath)
    easy_run.fully_buffered(cmd).raise_if_errors()
    assert os.path.exists('reference.json')
    ref_db = load.datablock('reference.json')[0]
    ref_imset = ref_db.extract_imagesets()[0]
    o = matrix.col(ref_imset.get_detector()[0].get_origin())
    ref_offsets = []
    for p in ref_imset.get_detector():
      intra_pnl = o - matrix.col(p.get_origin())
      ref_offsets.append(intra_pnl.length())
    assert approx_equal(offsets, ref_offsets)

  def tst_from_image_files(self):
    from subprocess import call, PIPE
    from glob import glob
    from os.path import join, exists

    # Find the image files
    image_files = glob(join(self.path, "centroid*.cbf"))
    image_files = ' '.join(image_files)

    # Import from the image files
    call('dials.import %s output.datablock=import_datablock.json' % \
           image_files, shell=True, stdout=PIPE)

    assert(exists("import_datablock.json"))


  def tst_from_template(self):
    from subprocess import call, PIPE
    from glob import glob
    from os.path import join, exists

    # Find the image files
    template = join(self.path, "centroid_####.cbf")

    # Import from the image files
    call('dials.import template=%s output.datablock=import_datablock.json' % \
           template, shell=True, stdout=PIPE)

    assert(exists("import_datablock.json"))


  def tst_extrapolate_scan(self):
    from subprocess import call, PIPE
    from os.path import join, exists

    # First image file
    image = join(self.path, "centroid_0001.cbf")

    cmd = 'dials.import %s output.datablock=import_extrapolate.json geometry.scan.image_range=1,900 geometry.scan.extrapolate_scan=True' % image

    # Import from the image file
    call(cmd, shell=True, stdout=PIPE)
    try:
      assert(exists("import_extrapolate.json"))
    except Exception:
      print cmd
      raise


if __name__ == '__main__':
  from dials.test import cd_auto
  with cd_auto(__file__):
    test = Test()
    test.run()

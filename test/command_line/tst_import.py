from __future__ import division
from libtbx.test_utils import approx_equal

class Test(object):
  def __init__(self):
    from os.path import join
    import libtbx.load_env
    try:
      dials_regression = libtbx.env.dist_path('dials_regression')
    except KeyError, e:
      print 'FAIL: dials_regression not configured'
      exit(0)

    self.path = join(dials_regression, "centroid_test_data")

  def run(self):
    self.tst_from_image_files()
    self.tst_import_beam_centre()
    self.tst_with_mask()
    self.tst_override_geometry()
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

    print 'OK'

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
    easy_run.fully_buffered(cmd)
    assert os.path.exists("datablock_with_mask.json")
    datablock = load.datablock("datablock_with_mask.json")[0]
    imgset = datablock.extract_imagesets()[0]
    assert imgset.external_lookup.mask.filename == mask_filename

    print 'OK'

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
            rotation_axis = 0,0,-1
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
    easy_run.fully_buffered(cmd)
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
    assert goniometer.get_rotation_axis() == (0,0,-1)
    assert goniometer.get_fixed_rotation() == (0,1,2,3,4,5,6,7,8)
    assert goniometer.get_setting_rotation() == (8,7,6,5,4,3,2,1,0)
    assert scan.get_image_range() == (1,4)
    assert scan.get_oscillation() == (1,2)

    print 'OK'

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
    easy_run.fully_buffered(cmd)
    assert os.path.exists("mosflm_beam_centre.json")
    datablock = load.datablock("mosflm_beam_centre.json")[0]
    imgset = datablock.extract_imagesets()[0]
    beam_centre = imgset.get_detector()[0].get_beam_centre(imgset.get_beam().get_s0())
    assert approx_equal(beam_centre, (200,100))

    # provide an alternative datablock.json to get geometry from
    cmd = 'dials.import %s reference_geometry=mosflm_beam_centre.json output.datablock=mosflm_beam_centre2.json' %image_files
    easy_run.fully_buffered(cmd)
    assert os.path.exists("mosflm_beam_centre2.json")
    datablock = load.datablock("mosflm_beam_centre2.json")[0]
    imgset = datablock.extract_imagesets()[0]
    beam_centre = imgset.get_detector()[0].get_beam_centre(imgset.get_beam().get_s0())
    assert approx_equal(beam_centre, (200,100))

    print 'OK'

  def tst_from_image_files(self):
    from subprocess import call, PIPE
    from glob import glob
    from os.path import join, exists

    # Find the image files
    image_files = glob(join(self.path, "centroid*.cbf"))
    image_files = ' '.join(image_files)

    # Import from the image files
    call('dials.import %s output.datablock=import_datablock.json' % image_files, shell=True, stdout=PIPE)

    # Get the expected output
    #expected = self.expected_import_from_image_files()

    assert(exists("import_datablock.json"))
    ## Read the created file and do a diff
    #with open("datablock.json", "r") as infile:
      #lines_a = infile.read().splitlines()
      #lines_a = [l.strip() for l in lines_a if "\"template\"" not in l]
      #diff = list(difflib.context_diff(
        #lines_a,
        #[l.strip() for l in expected.splitlines()]))
      #n = len(diff)
      #for i, line in enumerate(diff):
        #print line
      #assert(n == 0)

    print 'OK'

  #def expected_import_from_image_files(self):
    #return '''[
  #{
    #"__id__": "DataBlock",
    #"imageset": [
      #{
        #"__id__": "ImageSweep",
        #"beam": 0,
        #"detector": 0,
        #"goniometer": 0,
        #"scan": 0
      #}
    #],
    #"beam": [
      #{
        #"direction": [
          #0.0,
          #0.0,
          #1.0
        #],
        #"polarization_normal": [
          #0.0,
          #1.0,
          #0.0
        #],
        #"divergence": 0.0,
        #"polarization_fraction": 0.999,
        #"sigma_divergence": 0.0,
        #"wavelength": 0.9795
      #}
    #],
    #"detector": [
      #{
        #"panels": [
          #{
            #"origin": [
              #-212.47848,
              #220.00176,
              #-190.17999999999998
            #],
            #"fast_axis": [
              #1.0,
              #0.0,
              #0.0
            #],
            #"name": "Panel",
            #"slow_axis": [
              #0.0,
              #-1.0,
              #0.0
            #],
            #"mask": [
              #[
                #488,
                #1,
                #494,
                #2527
              #],
              #[
                #982,
                #1,
                #988,
                #2527
              #],
              #[
                #1476,
                #1,
                #1482,
                #2527
              #],
              #[
                #1970,
                #1,
                #1976,
                #2527
              #],
              #[
                #1,
                #196,
                #2463,
                #212
              #],
              #[
                #1,
                #408,
                #2463,
                #424
              #],
              #[
                #1,
                #620,
                #2463,
                #636
              #],
              #[
                #1,
                #832,
                #2463,
                #848
              #],
              #[
                #1,
                #1044,
                #2463,
                #1060
              #],
              #[
                #1,
                #1256,
                #2463,
                #1272
              #],
              #[
                #1,
                #1468,
                #2463,
                #1484
              #],
              #[
                #1,
                #1680,
                #2463,
                #1696
              #],
              #[
                #1,
                #1892,
                #2463,
                #1908
              #],
              #[
                #1,
                #2104,
                #2463,
                #2120
              #],
              #[
                #1,
                #2316,
                #2463,
                #2332
              #]
            #],
            #"trusted_range": [
              #-1.0,
              #495976.0
            #],
            #"image_size": [
              #2463,
              #2527
            #],
            #"px_mm_strategy": {
              #"type": "ParallaxCorrectedPxMmStrategy"
            #},
            #"type": "SENSOR_PAD",
            #"pixel_size": [
              #0.17200000000000001,
              #0.17200000000000001
            #]
          #}
        #]
      #}
    #],
    #"goniometer": [
      #{
        #"fixed_rotation": [
          #1.0,
          #0.0,
          #0.0,
          #0.0,
          #1.0,
          #0.0,
          #0.0,
          #0.0,
          #1.0
        #],
        #"rotation_axis": [
          #1.0,
          #0.0,
          #0.0
        #]
      #}
    #],
    #"scan": [
      #{
        #"exposure_time": [
          #0.2,
          #0.2,
          #0.2,
          #0.2,
          #0.2,
          #0.2,
          #0.2,
          #0.2,
          #0.2
        #],
        #"epochs": [
          #1360324992.0,
          #1360324992.0,
          #1360324993.0,
          #1360324993.0,
          #1360324993.0,
          #1360324993.0,
          #1360324993.0,
          #1360324994.0,
          #1360324994.0
        #],
        #"image_range": [
          #1,
          #9
        #],
        #"oscillation": [
          #0.0,
          #0.2
        #]
      #}
    #]
  #}
#]'''

  #def expected_import_from_xds_files(self):
    #return '''{
  #"__id__": "ExperimentList",
  #"experiment": [
    #{
      #"__id__": "Experiment",
      #"beam": 0,
      #"detector": 0,
      #"goniometer": 0,
      #"scan": 0,
      #"crystal": 0,
      #"imageset": 0
    #}
  #],
  #"imageset": [
    #{
      #"__id__": "ImageSweep",
    #}
  #],
  #"beam": [
    #{
      #"direction": [
        #-0.007852057721998333,
        #3.772524827250213e-14,
        #0.9999691721195861
      #],
      #"polarization_normal": [
        #0.0,
        #1.0,
        #0.0
      #],
      #"divergence": 0.0,
      #"polarization_fraction": 0.999,
      #"sigma_divergence": 0.058,
      #"wavelength": 0.9795
    #}
  #],
  #"detector": [
    #{
      #"panels": [
        #{
          #"origin": [
            #-211.53596470096178,
            #219.45303890619488,
            #-192.7062494437063
          #],
          #"fast_axis": [
            #0.9999551354884303,
            #0.0021159302715049923,
            #0.009233084500921031
          #],
          #"name": "Panel",
          #"slow_axis": [
            #0.0021250002879257116,
            #-0.999997269169901,
            #-0.0009726389448611214
          #],
          #"mask": [],
          #"trusted_range": [
            #0.0,
            #0.0
          #],
          #"image_size": [
            #2463,
            #2527
          #],
          #"px_mm_strategy": {
            #"type": "ParallaxCorrectedPxMmStrategy"
          #},
          #"type": "SENSOR_UNKNOWN",
          #"pixel_size": [
            #0.172,
            #0.172
          #]
        #}
      #]
    #}
  #],
  #"goniometer": [
    #{
      #"fixed_rotation": [
        #1.0,
        #0.0,
        #0.0,
        #0.0,
        #1.0,
        #0.0,
        #0.0,
        #0.0,
        #1.0
      #],
      #"rotation_axis": [
        #1.0,
        #-1.5919306617286774e-16,
        #-6.904199434387693e-16
      #]
    #}
  #],
  #"scan": [
    #{
      #"exposure_time": [
        #0.2,
        #0.2,
        #0.2,
        #0.2,
        #0.2,
        #0.2,
        #0.2,
        #0.2,
        #0.2
      #],
      #"epochs": [
        #1360324992.0,
        #1360324992.0,
        #1360324993.0,
        #1360324993.0,
        #1360324993.0,
        #1360324993.0,
        #1360324993.0,
        #1360324994.0,
        #1360324994.0
      #],
      #"image_range": [
        #1,
        #9
      #],
      #"oscillation": [
        #0.0,
        #0.2
      #]
    #}
  #],
  #"crystal": [
    #{
      #"__id__": "crystal",
      #"real_space_a": [
        #35.23781811553089,
        #-7.600614003857873,
        #22.077690418635804
      #],
      #"real_space_b": [
        #-22.657129890916668,
        #-1.4698317405529955,
        #35.65693038892429
      #],
      #"real_space_c": [
        #-5.295803077552249,
        #-38.99952334925477,
        #-4.972795822746061
      #],
      #"space_group_hall_symbol": " P 4 2",
      #"mosaicity": 0.157
    #}
  #]
#}'''

if __name__ == '__main__':
  from dials.test import cd_auto
  with cd_auto(__file__):
    test = Test()
    test.run()

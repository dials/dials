from __future__ import absolute_import, division, print_function

from glob import glob
import os

from procrunner import run_process
import pytest

def test_multiple_sweep_import_fails_without_allow_parameter(dials_regression, tmpdir):
  tmpdir.chdir()

  # Find the image files
  image_files = sorted(glob(os.path.join(dials_regression, "centroid_test_data", "centroid*.cbf")))
  del image_files[4] # Delete filename to force two sweeps

  # run without allowing multiple sweeps
  result = run_process(['dials.import'] + image_files + ['output.datablock=datablock_multiple_sweeps.json'])
  assert result['exitcode'] == 1
  assert 'ore than 1 sweep' in result['stderr']
  assert not os.path.exists("datablock_multiple_sweeps.json")

def test_multiple_sweep_import_suceeds_with_allow_parameter(dials_regression, tmpdir):
  tmpdir.chdir()

  # Find the image files
  image_files = sorted(glob(os.path.join(dials_regression, "centroid_test_data", "centroid*.cbf")))
  del image_files[4] # Delete filename to force two sweeps

  result = run_process(['dials.import'] + image_files + ['output.datablock=datablock_multiple_sweeps.json', 'allow_multiple_sweeps=True'])
  assert result['exitcode'] == 0
  assert result['stderr'] == ''
  assert os.path.exists("datablock_multiple_sweeps.json")

  from dxtbx.serialize import load
  datablock = load.datablock("datablock_multiple_sweeps.json")[0]
  imgset = datablock.extract_imagesets()
  assert len(imgset) == 2

def test_with_mask(dials_regression, tmpdir):
  tmpdir.chdir()

  # Find the image files
  image_files = glob(os.path.join(dials_regression, "centroid_test_data", "centroid*.cbf"))
  mask_filename = os.path.join(dials_regression, "centroid_test_data", "mask.pickle")

  result = run_process(['dials.import'] + image_files + ['mask=' + mask_filename, 'output.datablock=datablock_with_mask.json'])
  assert result['exitcode'] == 0
  assert result['stderr'] == ''
  assert os.path.exists("datablock_with_mask.json")

  from dxtbx.serialize import load
  datablock = load.datablock("datablock_with_mask.json")[0]
  imgset = datablock.extract_imagesets()[0]
  assert imgset.external_lookup.mask.filename == mask_filename

def test_override_geometry(dials_regression, tmpdir):
  tmpdir.chdir()

  # Find the image files
  image_files = glob(os.path.join(dials_regression, "centroid_test_data", "centroid*.cbf"))

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

  result = run_process(['dials.import'] + image_files + ['geometry.phil', 'output.datablock=override_geometry.json'])
  assert result['exitcode'] == 0
  assert result['stderr'] == ''
  assert os.path.exists("override_geometry.json")

  from dxtbx.serialize import load
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

def tst_import_beam_centre(dials_regression, tmpdir):
  tmpdir.chdir()

  # Find the image files
  image_files = glob(os.path.join(dials_regression, "centroid_test_data", "centroid*.cbf"))
  image_files = ' '.join(image_files)

  # provide mosflm beam centre to dials.import
  result = run_process(['dials.import', 'mosflm_beam_centre=100,200', 'output.datablock=mosflm_beam_centre.json'] + image_files)
  assert result['exitcode'] == 0
  assert result['stderr'] == ''
  assert os.path.exists("mosflm_beam_centre.json")

  from dxtbx.serialize import load
  datablock = load.datablock("mosflm_beam_centre.json")[0]
  imgset = datablock.extract_imagesets()[0]
  beam_centre = imgset.get_detector()[0].get_beam_centre(imgset.get_beam().get_s0())
  assert beam_centre == pytest.approx((200,100))

  # provide an alternative datablock.json to get geometry from
  result = run_process(['dials.import', 'reference_geometry=mosflm_beam_centre.json', 'output.datablock=mosflm_beam_centre2.json'] + image_files)
  assert result['exitcode'] == 0
  assert result['stderr'] == ''
  assert os.path.exists("mosflm_beam_centre2.json")
  datablock = load.datablock("mosflm_beam_centre2.json")[0]
  imgset = datablock.extract_imagesets()[0]
  beam_centre = imgset.get_detector()[0].get_beam_centre(imgset.get_beam().get_s0())
  assert beam_centre == pytest.approx((200,100))

def test_slow_fast_beam_centre(dials_regression, tmpdir):
  tmpdir.chdir()

  # test slow_fast_beam_centre with a multi-panel CS-PAD image
  impath = os.path.join(dials_regression, "image_examples",
      "LCLS_cspad_nexus", "idx-20130301060858401.cbf")
  result = run_process(['dials.import', 'slow_fast_beam_centre=134,42,18', 'output.datablock=slow_fast_beam_centre.json', impath])
  assert result['exitcode'] == 0
  assert result['stderr'] == ''
  assert os.path.exists('slow_fast_beam_centre.json')

  from dxtbx.serialize import load
  datablock = load.datablock('slow_fast_beam_centre.json')[0]
  imgset = datablock.extract_imagesets()[0]
  # beam centre on 18th panel
  s0 = imgset.get_beam().get_s0()
  beam_centre = imgset.get_detector()[18].get_beam_centre_px(s0)
  assert beam_centre == pytest.approx((42,134))

  # check relative panel positions have not changed
  from scitbx import matrix
  o = matrix.col(imgset.get_detector()[0].get_origin())
  offsets = []
  for p in imgset.get_detector():
    intra_pnl = o - matrix.col(p.get_origin())
    offsets.append(intra_pnl.length())

  result = run_process(['dials.import', 'output.datablock=reference.json', impath])
  assert result['exitcode'] == 0
  assert result['stderr'] == ''
  assert os.path.exists('reference.json')

  ref_db = load.datablock('reference.json')[0]
  ref_imset = ref_db.extract_imagesets()[0]
  o = matrix.col(ref_imset.get_detector()[0].get_origin())
  ref_offsets = []
  for p in ref_imset.get_detector():
    intra_pnl = o - matrix.col(p.get_origin())
    ref_offsets.append(intra_pnl.length())
  assert offsets == pytest.approx(ref_offsets)

def test_from_image_files(dials_regression, tmpdir):
  tmpdir.chdir()

  # Find the image files
  image_files = glob(os.path.join(dials_regression, "centroid_test_data", "centroid*.cbf"))

  # Import from the image files
  result = run_process(['dials.import'] + image_files + ['output.datablock=import_datablock.json'])
  assert result['exitcode'] == 0
  assert os.path.exists("import_datablock.json")

def test_from_template(dials_regression, tmpdir):
  tmpdir.chdir()

  # Find the image files
  template = os.path.join(dials_regression, "centroid_test_data", "centroid_####.cbf")

  # Import from the image files
  result = run_process(['dials.import', 'template=' + template, 'output.datablock=import_datablock.json'])
  assert result['exitcode'] == 0
  assert os.path.exists("import_datablock.json")

def test_extrapolate_scan(dials_regression, tmpdir):
  tmpdir.chdir()

  # First image file
  image = os.path.join(dials_regression, "centroid_test_data", "centroid_0001.cbf")

  result = run_process(['dials.import', image, 'output.datablock=import_extrapolate.json', 'geometry.scan.image_range=1,900', 'geometry.scan.extrapolate_scan=True'])
  assert result['exitcode'] == 0
  assert os.path.exists("import_extrapolate.json")

from __future__ import absolute_import, division, print_function

import glob
import os

import procrunner

def test_dials_show(dials_regression):
  path = os.path.join(dials_regression, "experiment_test_data", "experiment_1.json")
  result = procrunner.run(['dials.show', path], environment_override={'DIALS_NOBANNER': '1'})
  assert not result['exitcode'] and not result['stderr']
  output = list(filter(None, (s.rstrip() for s in result['stdout'].split('\n'))))
  assert "\n".join(output[4:]) == """
Experiment 0:
Detector:
Panel:
  name: Panel
  type: SENSOR_PAD
  identifier:
  pixel_size:{0.172,0.172}
  image_size: {2463,2527}
  trusted_range: {-1,495976}
  thickness: 0
  material:
  mu: 0
  gain: 1
  pedestal: 0
  fast_axis: {1,0,0}
  slow_axis: {0,-1,0}
  origin: {-212.478,220.002,-190.18}
  distance: 190.18
  pixel to millimeter strategy: SimplePxMmStrategy
Max resolution (at corners): 1.008178
Max resolution (inscribed):  1.204283
Beam:
    wavelength: 0.9795
    sample to source direction : {0,0,1}
    divergence: 0
    sigma divergence: 0
    polarization normal: {0,1,0}
    polarization fraction: 0.999
Beam centre (mm): (212.48,220.00)
Beam centre (px): (1235.34,1279.08)
Scan:
    image range:   {1,9}
    oscillation:   {0,0.2}
    exposure time: 0.2
Goniometer:
    Rotation axis:   {1,0,0}
    Fixed rotation:  {1,0,0,0,1,0,0,0,1}
    Setting rotation:{1,0,0,0,1,0,0,0,1}
Crystal:
    Unit cell: (42.272, 42.272, 39.670, 90.000, 89.999, 90.000)
    Space group: P 4 2 2
    U matrix:  {{ 0.8336, -0.5360, -0.1335},
                {-0.1798, -0.0348, -0.9831},
                { 0.5223,  0.8435, -0.1254}}
    B matrix:  {{ 0.0237,  0.0000,  0.0000},
                {-0.0000,  0.0237,  0.0000},
                {-0.0000,  0.0000,  0.0252}}
    A = UB:    {{ 0.0197, -0.0127, -0.0034},
                {-0.0043, -0.0008, -0.0248},
                { 0.0124,  0.0200, -0.0032}}
    Mosaicity:  0.157000
""".strip()

def test_dials_show_i04_weak_data(dials_regression):
  path = os.path.join(
    dials_regression, "indexing_test_data", "i04_weak_data", "datablock_orig.json")
  result = procrunner.run(["dials.show", path], environment_override={'DIALS_NOBANNER': '1'})
  assert not result['exitcode'] and not result['stderr']
  output = list(filter(None, (s.rstrip() for s in result['stdout'].split('\n'))))
  assert "\n".join(output[4:]) == """
Experiment 0:
Detector:
Panel:
  name: Panel
  type: SENSOR_PAD
  identifier:
  pixel_size:{0.172,0.172}
  image_size: {2463,2527}
  trusted_range: {-1,161977}
  thickness: 0
  material:
  mu: 0
  gain: 1
  pedestal: 0
  fast_axis: {1,0,0}
  slow_axis: {0,-1,0}
  origin: {-210.76,205.277,-265.27}
  distance: 265.27
  pixel to millimeter strategy: SimplePxMmStrategy
Max resolution (at corners): 1.161261
Max resolution (inscribed):  1.509475
Beam:
    wavelength: 0.97625
    sample to source direction : {0,0,1}
    divergence: 0
    sigma divergence: 0
    polarization normal: {0,1,0}
    polarization fraction: 0.999
Beam centre (mm): (210.76,205.28)
Beam centre (px): (1225.35,1193.47)
Scan:
    image range:   {1,540}
    oscillation:   {82,0.15}
    exposure time: 0.067
Goniometer:
    Rotation axis:   {1,0,0}
    Fixed rotation:  {1,0,0,0,1,0,0,0,1}
    Setting rotation:{1,0,0,0,1,0,0,0,1}
""".strip()

def test_dials_show_centroid_test_data(dials_regression):
  path = os.path.join(
    dials_regression, "centroid_test_data", "centroid_*.cbf")
  g = glob.glob(path)
  assert g, path
  result = procrunner.run(["dials.show"] + g, environment_override={'DIALS_NOBANNER': '1'})
  assert not result['exitcode'] and not result['stderr']
  output = list(filter(None, (s.rstrip() for s in result['stdout'].split('\n'))))
  assert "\n".join(output[4:]) == """
Experiment 0:
Detector:
Panel:
  name: Panel
  type: SENSOR_PAD
  identifier:
  pixel_size:{0.172,0.172}
  image_size: {2463,2527}
  trusted_range: {-1,495976}
  thickness: 0.32
  material: Si
  mu: 3.96038
  gain: 1
  pedestal: 0
  fast_axis: {1,0,0}
  slow_axis: {0,-1,0}
  origin: {-212.478,220.002,-190.18}
  distance: 190.18
  pixel to millimeter strategy: ParallaxCorrectedPxMmStrategy
    mu: 3.96038
    t0: 0.32
Max resolution (at corners): 1.008375
Max resolution (inscribed):  1.204621
Beam:
    wavelength: 0.9795
    sample to source direction : {0,0,1}
    divergence: 0
    sigma divergence: 0
    polarization normal: {0,1,0}
    polarization fraction: 0.999
Beam centre (mm): (212.48,220.00)
Beam centre (px): (1235.34,1279.08)
Scan:
    image range:   {1,9}
    oscillation:   {0,0.2}
    exposure time: 0.2
Goniometer:
    Rotation axis:   {1,0,0}
    Fixed rotation:  {1,0,0,0,1,0,0,0,1}
    Setting rotation:{1,0,0,0,1,0,0,0,1}
""".strip()

def test_dials_show_reflection_table(dials_regression):
  """Test the output of dials.show on a reflection_table pickle file"""
  path = os.path.join(dials_regression, "centroid_test_data", "integrated.pickle")
  result = procrunner.run(["dials.show", path], environment_override={'DIALS_NOBANNER': '1'})
  assert not result['exitcode'] and not result['stderr']
  output = list(filter(None, (s.rstrip() for s in result['stdout'].split('\n'))))
  assert output[4] == 'Reflection list contains 2269 reflections'
  headers = ['Column', 'min', 'max', 'mean']
  for header in headers:
    assert header in output[6]
  row_names = ['background.mean', 'background.sum.value',
    'background.sum.variance', 'd', 'dqe', 'flags', 'id', 'imageset_id',
    'intensity.prf.value', 'intensity.prf.variance', 'intensity.sum.value',
    'intensity.sum.variance', 'lp', 'miller_index', 'num_pixels.background',
    'num_pixels.background_used', 'num_pixels.foreground', 'num_pixels.valid',
    'panel', 'partial_id', 'partiality', 'profile.correlation', 'profile.rmsd',
    'rlp', 's1', 'shoebox', 'summed I', 'N pix', 'N valid foreground pix',
    'xyzcal.mm', 'xyzcal.px', 'xyzobs.mm.value', 'xyzobs.mm.variance',
    'xyzobs.px.value', 'xyzobs.px.variance', 'zeta']
  for (name, out) in zip(row_names, output[8:-1]):
    assert name in out

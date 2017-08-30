from __future__ import absolute_import, division

def run():
  import os
  import libtbx.load_env
  from libtbx import easy_run
  from libtbx.test_utils import show_diff
  try:
    dials_regression = libtbx.env.dist_path('dials_regression')
  except KeyError, e:
    print 'FAIL: dials_regression not configured'
    exit(0)

  path = os.path.join(dials_regression, "experiment_test_data")

  cmd = "dials.show %s/experiment_1.json" %path
  result = easy_run.fully_buffered(cmd).raise_if_errors()
  assert not show_diff("\n".join(result.stdout_lines[6:]), """\
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
  fast_axis: {1,0,0}
  slow_axis: {0,-1,0}
  origin: {-212.478,220.002,-190.18}
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
""", strip_trailing_whitespace=True)

  path = os.path.join(
    dials_regression, "indexing_test_data", "i04_weak_data")
  cmd = "dials.show %s/datablock_orig.json" %path
  result = easy_run.fully_buffered(cmd).raise_if_errors()
  assert not show_diff("\n".join(result.stdout_lines[8:]), """\
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
  fast_axis: {1,0,0}
  slow_axis: {0,-1,0}
  origin: {-210.76,205.277,-265.27}
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
""", strip_trailing_whitespace=True)

  path = os.path.join(
    dials_regression, "centroid_test_data", "centroid_*.cbf")
  import glob
  g = glob.glob(path)
  assert len(g) > 0, path
  cmd = "dials.show %s" %(' '.join(g))
  result = easy_run.fully_buffered(cmd).raise_if_errors()
  assert (
    "Format: <class 'dxtbx.format.FormatCBFMiniPilatus.FormatCBFMiniPilatus'>"
    in result.stdout_lines), result.show_stdout()
  assert not show_diff("\n".join(result.stdout_lines[8:]), """\
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
  fast_axis: {1,0,0}
  slow_axis: {0,-1,0}
  origin: {-212.478,220.002,-190.18}
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
""", strip_trailing_whitespace=True)


if __name__ == '__main__':
  from dials.test import cd_auto
  with cd_auto(__file__):
    run()
    print "OK"

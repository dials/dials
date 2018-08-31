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
  assert "\n".join(output[6:]) == """
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
  assert (
    "Format: <class 'dxtbx.format.FormatCBFMiniPilatus.FormatCBFMiniPilatus'>"
    in result['stdout'])
  output = list(filter(None, (s.rstrip() for s in result['stdout'].split('\n'))))
  assert "\n".join(output[6:]) == """
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
  assert "\n".join(output[4:]) == """
Reflection list contains 2269 reflections
---------------------------------------------------------------------------------------------------------------------------------------------
| Column                     | min                                | max                                | mean                               |
---------------------------------------------------------------------------------------------------------------------------------------------
| background.mean            | 0.0                                | 3.6                                | 0.5                                |
| background.sum.value       | 0.0                                | 925.1                              | 226.1                              |
| background.sum.variance    | 0.0                                | 1143.6                             | 265.7                              |
| d                          | 0.00                               | 13.68                              | 1.50                               |
| dqe                        | 0.000                              | 0.904                              | 0.709                              |
| flags                      | 769                                | 2228260                            | 502517                             |
| id                         | -1                                 | 0                                  | 0                                  |
| imageset_id                | 0                                  | 0                                  | 0                                  |
| intensity.prf.value        | -11.3                              | 40839.8                            | 426.6                              |
| intensity.prf.variance     | -1.0                               | 41045.5                            | 518.2                              |
| intensity.sum.value        | -52.4                              | 41476.6                            | 547.4                              |
| intensity.sum.variance     | 0.0                                | 42223.9                            | 816.1                              |
| lp                         | 0.000                              | 0.920                              | 0.486                              |
| miller_index               | -31, -36, -28                      | 6, 20, 32                          | -7, -4, 1                          |
| num_pixels.background      | 0                                  | 15219                              | 3195                               |
| num_pixels.background_used | 0                                  | 15219                              | 3195                               |
| num_pixels.foreground      | 0                                  | 1809                               | 484                                |
| num_pixels.valid           | 0                                  | 17028                              | 3680                               |
| panel                      | 0                                  | 0                                  | 0                                  |
| partial_id                 | 0                                  | 1990                               | 873                                |
| partiality                 | 0.0000                             | 0.9998                             | 0.7450                             |
| profile.correlation        | -0.235                             | 0.912                              | 0.199                              |
| profile.rmsd               | 0.000                              | 0.228                              | 0.101                              |
| rlp                        | -0.6434, -0.5749, 0.0000           | 0.5587, 0.6052, 0.2326             | -0.0042, 0.0034, 0.0121            |
| s1                         | -0.7430, -0.7560, -1.0209          | 0.7523, 0.7665, -0.5510            | -0.0042, 0.0055, -0.8036           |
| shoebox                    |                                    |                                    |                                    |
|   summed I                 | 0.0                                | 10514.0                            | 73.5                               |
|   N pix                    | 0.0                                | 378.0                              | 2.4                                |
|   N valid foreground pix   | 0.0                                | 87.0                               | 1.0                                |
| xyzcal.mm                  | 0.00, 0.00, -0.04                  | 422.90, 434.48, 0.08               | 204.72, 210.72, 0.01               |
| xyzcal.px                  | -0.23, 0.00, -11.21                | 2459.34, 2526.64, 23.88            | 1190.23, 1225.10, 4.18             |
| xyzobs.mm.value            | 0.61, 0.79, 0.00                   | 422.67, 433.87, 0.03               | 212.02, 218.32, 0.02               |
| xyzobs.mm.variance         | 0.0000e+00, 0.0000e+00, 0.0000e+00 | 7.9419e-02, 1.0688e-02, 1.5462e-05 | 3.4148e-03, 3.3928e-03, 1.2491e-06 |
| xyzobs.px.value            | 2.94, 3.96, 0.50                   | 2458.01, 2523.09, 8.50             | 1232.67, 1269.32, 4.46             |
| xyzobs.px.variance         | 0.0000, 0.0000, 0.0000             | 2.6845, 0.3613, 1.2690             | 0.1154, 0.1147, 0.1025             |
| zeta                       | -1.000                             | 1.000                              | -0.003                             |
---------------------------------------------------------------------------------------------------------------------------------------------""".strip()

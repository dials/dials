from __future__ import absolute_import, division, print_function

from math import pi

def test_polarization_conversion():
  from random import uniform
  from scitbx import matrix
  from dials.util.nexus.nx_mx import polarization_normal_to_stokes
  from dials.util.nexus.nx_mx import polarization_stokes_to_normal

  EPS = 1e-7

  def conv(n, p):
    S0, S1, S2, S3 = polarization_normal_to_stokes(n, p)
    return polarization_stokes_to_normal(S0, S1, S2, S3)

  for i in range(1000):

    n1 = matrix.col((uniform(-1,1), uniform(-1,1), 0)).normalize()
    p1 = uniform(0,1)
    n2, p2 = conv(n1, p1)
    angle = n2.angle(n1)
    assert(abs(angle) < EPS or abs(angle - pi) < EPS or abs(angle-2*pi) < EPS)
    assert(abs(p1 - p2) < EPS)

def run_single(experiments1, filename):
  from dials.util.nexus import dump, load
  from scitbx import matrix

  EPS = 1e-7

  try:
    run_single.count += 1
  except Exception:
    run_single.count = 0

  # Dump the file
  dump(experiments1, None, filename)

  # Load the file
  experiments2, reflections = load(filename)
  assert(experiments2 is not None)
  assert(reflections is None)
  assert(len(experiments2) == len(experiments1))

  index1 = []
  index2 = []

  for exp1, exp2 in zip(experiments1, experiments2):

    # Check the beam
    b1 = exp1.beam
    b2 = exp2.beam
    assert(all(abs(d1 - d2) < EPS
      for d1, d2 in zip(b1.get_direction(), b2.get_direction())))
    assert(abs(b1.get_wavelength() - b2.get_wavelength()) < EPS)
    assert(abs(b1.get_polarization_fraction() - b2.get_polarization_fraction()) < EPS)
    n1 = matrix.col(b1.get_polarization_normal())
    n2 = matrix.col(b2.get_polarization_normal())
    angle = n2.angle(n1)
    assert(abs(angle) < EPS or abs(angle - pi) < EPS or abs(angle-2*pi) < EPS)

    # Check the goniometer
    g1 = exp1.goniometer
    g2 = exp2.goniometer
    if g1 is not None:
      assert(all(abs(d1 - d2) < EPS
        for d1, d2 in zip(g1.get_rotation_axis(), g2.get_rotation_axis())))
      assert(all(abs(d1 - d2) < EPS
        for d1, d2 in zip(g1.get_fixed_rotation(), g2.get_fixed_rotation())))
      assert(all(abs(d1 - d2) < EPS
        for d1, d2 in zip(g1.get_setting_rotation(), g2.get_setting_rotation())))
    else:
      assert(g2 is None)

    # Check the scan
    s1 = exp1.scan
    s2 = exp2.scan
    if s1 is not None:
      assert(len(s1) == len(s2))
      assert(s1.get_image_range() == s2.get_image_range())
      assert(abs(s1.get_oscillation()[0] - s1.get_oscillation()[0]) < EPS)
      assert(abs(s1.get_oscillation()[1] - s1.get_oscillation()[1]) < EPS)
      for e1, e2 in zip(s1.get_exposure_times(), s2.get_exposure_times()):
        assert(abs(e1 - e2) < EPS)
      for e1, e2 in zip(s1.get_epochs(), s2.get_epochs()):
        assert(abs(e1 - e2) < EPS)
    else:
      assert(s2 is None)

    # Check the detector
    d1 = exp1.detector
    d2 = exp2.detector
    assert(len(d1) == len(d2))
    for p1, p2 in zip(d1, d2):
      assert(p1.get_type() == p2.get_type())
      assert(p1.get_material() == p2.get_material())
      assert(p1.get_thickness() == p2.get_thickness())
      assert(p1.get_image_size() == p2.get_image_size())
      assert(p1.get_pixel_size() == p2.get_pixel_size())
      assert(p1.get_trusted_range() == p2.get_trusted_range())
      for x1, x2 in zip(p1.get_fast_axis(), p2.get_fast_axis()):
        assert(abs(x1 - x2) < EPS)
      for x1, x2 in zip(p1.get_slow_axis(), p2.get_slow_axis()):
        assert(abs(x1 - x2) < EPS)
      for x1, x2 in zip(p1.get_origin(), p2.get_origin()):
        assert(abs(x1 - x2) < EPS)

    # Check the crystal
    c1 = exp1.crystal
    c2 = exp2.crystal
    assert(c1.get_space_group() == c2.get_space_group())
    for p1, p2 in zip(c1.get_unit_cell().parameters(),
                      c2.get_unit_cell().parameters()):
      assert(abs(p1 -p2) < EPS)
    for p1, p2 in zip(c1.get_A(),
                      c2.get_A()):
      assert(abs(p1 -p2) < EPS)
    assert(c1.num_scan_points == c2.num_scan_points)
    for i in range(c1.num_scan_points):
      A1 = c1.get_A_at_scan_point(i)
      A2 = c2.get_A_at_scan_point(i)
      for a1, a2 in zip(A1, A2):
        assert(abs(a1 - a2) < EPS)
      uc1 = c1.get_unit_cell_at_scan_point(i)
      uc2 = c2.get_unit_cell_at_scan_point(i)
      for p1, p2 in zip(uc1.parameters(), uc2.parameters()):
        assert(abs(p1-p2) < EPS)

    index1.append((
      id(exp1.beam),
      id(exp1.detector),
      id(exp1.goniometer),
      id(exp1.scan),
      id(exp1.crystal)))

    index2.append((
      id(exp2.beam),
      id(exp2.detector),
      id(exp2.goniometer),
      id(exp2.scan),
      id(exp2.crystal)))

  # Get a list of all beam etc
  beam1, detector1, goniometer1, scan1, crystal1 = zip(*index1)
  beam2, detector2, goniometer2, scan2, crystal2 = zip(*index2)

  # If any models are shared then check they are shard in both
  num = len(beam1)
  for i in range(0,num-1):
    for j in range(1,num):
      if beam1[i] == beam1[j]:
        assert(beam2[i] == beam2[j])
      else:
        assert(beam2[i] != beam2[j])
      if detector1[i] == detector1[j]:
        assert(detector2[i] == detector2[j])
      else:
        assert(detector2[i] != detector2[j])
      if goniometer1[i] == goniometer1[j]:
        assert(goniometer2[i] == goniometer2[j])
      else:
        assert(goniometer2[i] != goniometer2[j])
      if scan1[i] == scan1[j]:
        assert(scan2[i] == scan2[j])
      else:
        assert(scan2[i] != scan2[j])
      if crystal1[i] == crystal1[j]:
        assert(crystal2[i] == crystal2[j])
      else:
        assert(crystal2[i] != crystal2[j])

def test_run(run_in_tmpdir, dials_regression):
  from dxtbx.model.experiment_list import ExperimentListFactory
  from os.path import join
  path = join(dials_regression, "nexus_test_data", "shared_models")
  filename_list = [
    'single',
    'multiple_unrelated',
    'multi_crystal',
    'two_colour',
    'multiple_sweeps',
    'stills'
  ]
  for filename in filename_list:
    filename_in = join(path, "%s.json" % filename)
    filename_out = "%s.nxs" % filename
    experiments = ExperimentListFactory.from_json_file(filename_in)
    run_single(experiments, filename_out)

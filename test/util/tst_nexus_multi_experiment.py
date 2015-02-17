
from __future__ import division

def run_single(experiments1):
  from dials.util.nexus import dump, load

  EPS = 1e-7

  try:
    run_single.count += 1
  except Exception:
    run_single.count = 0

  # Dump the file
  dump(experiments1, None, "hklout_%d.nxs" % run_single.count)

  # Load the file
  experiments2, reflections = load("hklout_%d.nxs" % run_single.count)
  assert(experiments2 is not None)
  assert(reflections is None)
  assert(len(experiments2) == len(experiments1))

  for exp1, exp2 in zip(experiments1, experiments2):

    # Check the beam
    b1 = exp1.beam
    b2 = exp2.beam
    assert(all(abs(d1 - d2) < EPS
      for d1, d2 in zip(b1.get_direction(), b2.get_direction())))
    assert(abs(b1.get_wavelength() - b2.get_wavelength()) < EPS)
    assert(abs(b1.get_polarization_fraction() - b2.get_polarization_fraction()) < EPS)
    assert(all(abs(d1 - d2) < EPS
      for d1, d2 in zip(b1.get_polarization_normal(), b2.get_polarization_normal())))
    print 'OK'

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
    print 'OK'

    # Check the scan
    s1 = exp1.scan
    s2 = exp2.scan
    print s1
    print s2
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
    print 'OK'

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
    print 'OK'

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
    print 'OK'

  print 'OK'

def run():
  from dxtbx.model.experiment.experiment_list import ExperimentListFactory
  from os.path import join
  import libtbx.load_env
  try:
    dials_regression = libtbx.env.dist_path('dials_regression')
  except KeyError, e:
    print 'FAIL: dials_regression not configured'
    exit(0)
  path = join(dials_regression, "nexus_test_data", "shared_models")
  filename_list = [
    'single.json',
    'multiple_unrelated.json',
    'multi_crystal.json',
    'two_colour.json',
    'multiple_sweeps.json',
    'stills.json'
  ]
  for filename in filename_list:
    experiments = ExperimentListFactory.from_json_file(join(path, filename))
    print filename
    run_single(experiments)

if __name__ == '__main__':
  from dials.test import cd_auto
  with cd_auto(__file__):
    run()

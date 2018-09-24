from __future__ import absolute_import, division, print_function

def test_run(dials_regression, run_in_tmpdir):
  from dials.util.nexus import dump, load
  from dxtbx.model.experiment_list import ExperimentListFactory
  from dials.array_family import flex
  from os.path import join
  path = join(dials_regression, "nexus_test_data")

  # Read the experiments
  experiments1 = ExperimentListFactory.from_json_file(
    join(path, "refined_experiments.json"))

  # Read the reflections
  reflections1 = flex.reflection_table.from_pickle(
    join(path, "integrated.pickle"))

  # Delete some columns for the test
  del reflections1['s1']
  del reflections1['zeta']
  del reflections1['background.mse']

  # Dump the reflections
  dump(experiments1, reflections1, "hklout.nxs")

  # Load them again
  experiments2, reflections2 = load("hklout.nxs")

  EPS = 1e-7

  # Check the reflections are OK
  assert(reflections1.nrows() == reflections2.nrows())
  assert(reflections1.ncols() == reflections2.ncols())
  for key in reflections1.keys():
    data1 = reflections1[key]
    data2 = reflections2[key]
    assert(data1.__class__ == data2.__class__)
    if isinstance(data1, flex.double):
      assert(data1.all_approx_equal(data2))
    elif isinstance(data1, flex.int6):
      for p1, p2 in zip(data1.parts(), data2.parts()):
        assert(p1.all_eq(p2))
    elif isinstance(data1, flex.vec3_double):
      for p1, p2 in zip(data1.parts(), data2.parts()):
        assert(p1.all_approx_equal(p2))
    else:
      assert(data1.all_eq(data2))

  # Test passed

  # Check the experiments are ok
  assert(len(experiments1) == len(experiments2))
  exp1 = experiments1[0]
  exp2 = experiments2[0]

  # Check the beam
  b1 = exp1.beam
  b2 = exp2.beam
  assert(all(abs(d1 - d2) < EPS
    for d1, d2 in zip(b1.get_direction(), b2.get_direction())))
  assert(abs(b1.get_wavelength() - b2.get_wavelength()) < EPS)
  assert(abs(b1.get_polarization_fraction() - b2.get_polarization_fraction()) < EPS)
  assert(all(abs(d1 - d2) < EPS
    for d1, d2 in zip(b1.get_polarization_normal(), b2.get_polarization_normal())))

  # Check the goniometer
  g1 = exp1.goniometer
  g2 = exp2.goniometer
  assert(all(abs(d1 - d2) < EPS
    for d1, d2 in zip(g1.get_rotation_axis(), g2.get_rotation_axis())))
  assert(all(abs(d1 - d2) < EPS
    for d1, d2 in zip(g1.get_fixed_rotation(), g2.get_fixed_rotation())))
  assert(all(abs(d1 - d2) < EPS
    for d1, d2 in zip(g1.get_setting_rotation(), g2.get_setting_rotation())))

  # Check the scan
  s1 = exp1.scan
  s2 = exp2.scan
  assert(len(s1) == len(s2))
  assert(s1.get_image_range() == s2.get_image_range())
  assert(abs(s1.get_oscillation()[0] - s1.get_oscillation()[0]) < EPS)
  assert(abs(s1.get_oscillation()[1] - s1.get_oscillation()[1]) < EPS)
  for e1, e2 in zip(s1.get_exposure_times(), s2.get_exposure_times()):
    assert(abs(e1 - e2) < EPS)
  for e1, e2 in zip(s1.get_epochs(), s2.get_epochs()):
    assert(abs(e1 - e2) < EPS)

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

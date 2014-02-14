from __future__ import division

class Test(object):

  def __init__(self):
    import os
    import libtbx.load_env
    from cctbx.crystal.crystal_model.serialize import load_crystal
    from dials.model.serialize import load
    from dials.algorithms.shoebox import BBoxCalculator
    from dials.algorithms.shoebox import MaskForeground
    from dials.model.experiment.experiment_list import Experiment

    try:
      dials_regression = libtbx.env.dist_path('dials_regression')
    except KeyError, e:
      print 'FAIL: dials_regression not configured'
      exit(0)

    # Set the sweep filename and load the sweep
    sweep_filename = os.path.join(dials_regression, 'centroid_test_data',
        'sweep.json')
    crystal_filename = os.path.join(dials_regression, 'centroid_test_data',
        'crystal.json')

    # Load the sweep
    self.sweep = load.sweep(sweep_filename)
    self.crystal = load_crystal(crystal_filename)
    self.beam = self.sweep.get_beam()
    self.detector = self.sweep.get_detector()
    self.goniometer = self.sweep.get_goniometer()
    self.scan = self.sweep.get_scan()
    self.experiment = Experiment(
      imageset = self.sweep,
      beam = self.beam,
      detector = self.detector,
      goniometer = self.goniometer,
      scan = self.scan,
      crystal = self.crystal)
    self.delta_d = 5 * self.beam.get_sigma_divergence(deg=False)
    self.delta_m = 5 * self.crystal.get_mosaicity(deg=False)
    self.nsigma = 5
    assert(len(self.detector) == 1)

    # Get the function object to mask the foreground
    self.mask_foreground = MaskForeground(self.beam, self.detector,
        self.goniometer, self.scan, self.delta_d, self.delta_m)


  def run(self):
    from scitbx import matrix
    from dials.algorithms.shoebox import MaskCode
    from dials.algorithms.reflection_basis import CoordinateSystem
    assert(len(self.detector) == 1)
    s0 = self.beam.get_s0()
    m2 = self.goniometer.get_rotation_axis()
    s0_length = matrix.col(self.beam.get_s0()).length()
    width, height = self.detector[0].get_image_size()

    # Generate some reflections
    reflections = self.generate_reflections(100)

    # Mask the foreground in each
    self.mask_foreground(
      reflections['shoebox'],
      reflections['s1'],
      reflections['xyzcal.px'].parts()[2])

    # Loop through all the reflections and check the mask values
    shoebox = reflections['shoebox']
    beam_vector = reflections['s1']
    rotation_angle = reflections['xyzcal.mm'].parts()[2]
    for i in range(len(reflections)):
      mask = shoebox[i].mask[0:1,:,:]
      x0, x1, y0, y1, z0, z1 = shoebox[i].bbox
      s1 = beam_vector[i]
      phi = rotation_angle[i]
      cs = CoordinateSystem(m2, s0, s1, phi)
      for j in range(y1 - y0):
        for i in range(x1 - x0):
          value1 = mask[0, j, i]
          s1 = self.detector[0].get_pixel_lab_coord(
              (x0 + i + 0.5, y0 + j + 0.5))
          s1 = matrix.col(s1).normalize() * s0_length
          e1, e2 = cs.from_beam_vector(s1)
          aa = (e1 / self.delta_d)**2
          bb = (e2 / self.delta_d)**2
          if (x0 + i < 0 or y0 + j < 0 or
              x0 + i > width or y0 + j > height):
            value2 = MaskCode.Valid
          else:
            if (aa + bb <= 1.0):
              value2 = MaskCode.Valid | MaskCode.Foreground
            else:
              value2 = MaskCode.Valid | MaskCode.Background
          assert(value1 == value2)

    # Test passed
    print 'OK'

  def generate_reflections(self, num):
    from random import randint
    from scitbx import matrix
    from dials.array_family import flex
    assert(len(self.detector) == 1)
    beam_vector = flex.vec3_double(num)
    xyzcal_px = flex.vec3_double(num)
    xyzcal_mm = flex.vec3_double(num)
    panel = flex.size_t(num)
    s0_length = matrix.col(self.beam.get_s0()).length()
    for i in range(num):
      x = randint(0, 2000)
      y = randint(0, 2000)
      z = randint(0, 8)
      s1 = self.detector[0].get_pixel_lab_coord((x, y))
      s1 = matrix.col(s1).normalize() * s0_length
      phi = self.scan.get_angle_from_array_index(z, deg=False)
      beam_vector[i] = s1
      xyzcal_px[i] = (x, y, z)
      (x, y)  = self.detector[0].pixel_to_millimeter((x, y))
      xyzcal_mm[i] = (x, y, phi)
      panel[i] = 0

    rlist = flex.reflection_table()
    rlist['s1'] = beam_vector
    rlist['panel'] = panel
    rlist['xyzcal.px'] = xyzcal_px
    rlist['xyzcal.mm'] = xyzcal_mm
    rlist.compute_bbox(self.experiment, self.nsigma)
    rlist['shoebox'] = flex.shoebox(
      rlist['panel'], rlist['bbox'])
    rlist['shoebox'].allocate()
    return rlist

if __name__ == '__main__':
  test = Test()
  test.run()

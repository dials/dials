from __future__ import division

class Test(object):

  def __init__(self):
    import os
    import libtbx.load_env
    from dxtbx.serialize.load import crystal as load_crystal
    from dials.model.serialize import load
    from dials.algorithms.profile_model.gaussian_rs import Model
    from dials.algorithms.profile_model.gaussian_rs import MaskCalculator3D
    from dxtbx.model.experiment.experiment_list import Experiment, ExperimentList

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
    self.delta_d = 3 * self.beam.get_sigma_divergence(deg=False)
    self.delta_m = 3 * self.crystal.get_mosaicity(deg=False)
    self.nsigma = 3
    self.profile_model = Model(
      None,
      self.nsigma,
      self.beam.get_sigma_divergence(deg=False),
      self.crystal.get_mosaicity(deg=False))
    self.experiment = ExperimentList()
    self.experiment.append(Experiment(
      imageset = self.sweep,
      beam = self.beam,
      detector = self.detector,
      goniometer = self.goniometer,
      scan = self.scan,
      crystal = self.crystal,
      profile = self.profile_model))


    assert(len(self.detector) == 1)

    # Get the function object to mask the foreground
    self.mask_foreground = MaskCalculator3D(self.beam, self.detector,
        self.goniometer, self.scan, self.delta_d, self.delta_m)


  def run(self):
    from scitbx import matrix
    from scitbx.array_family import flex
    from dials.algorithms.shoebox import MaskCode
    from dials.algorithms.profile_model.gaussian_rs import CoordinateSystem
    from math import sqrt
    assert(len(self.detector) == 1)
    s0 = self.beam.get_s0()
    m2 = self.goniometer.get_rotation_axis()
    s0_length = matrix.col(self.beam.get_s0()).length()
    width, height = self.detector[0].get_image_size()
    zrange = self.scan.get_array_range()
    phi0, dphi = self.scan.get_oscillation(deg=False)

    # Generate some reflections
    reflections = self.generate_reflections(10)

    # Mask the foreground in each
    self.mask_foreground(
      reflections['shoebox'],
      reflections['s1'],
      reflections['xyzcal.px'].parts()[2],
      reflections['panel'])

    # Loop through all the reflections and check the mask values
    shoebox = reflections['shoebox']
    beam_vector = reflections['s1']
    rotation_angle = reflections['xyzcal.mm'].parts()[2]
    for l in range(len(reflections)):
      mask = shoebox[l].mask
      x0, x1, y0, y1, z0, z1 = shoebox[l].bbox
      s1 = beam_vector[l]
      phi = rotation_angle[l]
      cs = CoordinateSystem(m2, s0, s1, phi)

      def rs_coord(i, j, k):
        s1d = self.detector[0].get_pixel_lab_coord((i,j))
        s1d = matrix.col(s1d).normalize() * s0_length
        e1, e2 = cs.from_beam_vector(s1d)
        e3 = cs.from_rotation_angle_fast(phi0 + (k - zrange[0]) * dphi)
        return e1, e2, e3

      new_mask = flex.int(mask.accessor(), 0)
      for k in range(z1 - z0):
        for j in range(y1 - y0):
          for i in range(x1 - x0):
            #value1 = mask[k, j, i]
            e11, e12, e13 = rs_coord(x0 + i,     y0 + j,     z0 + k)
            e21, e22, e23 = rs_coord(x0 + i + 1, y0 + j,     z0 + k)
            e31, e32, e33 = rs_coord(x0 + i,     y0 + j + 1, z0 + k)
            e41, e42, e43 = rs_coord(x0 + i,     y0 + j,     z0 + k + 1)
            e51, e52, e53 = rs_coord(x0 + i + 1, y0 + j + 1, z0 + k)
            e61, e62, e63 = rs_coord(x0 + i + 1, y0 + j,     z0 + k + 1)
            e71, e72, e73 = rs_coord(x0 + i,     y0 + j + 1, z0 + k + 1)
            e81, e82, e83 = rs_coord(x0 + i + 1, y0 + j + 1, z0 + k + 1)
            de1 = (e11/self.delta_d)**2+(e12/self.delta_d)**2#+(e13/self.delta_m)**2
            de2 = (e21/self.delta_d)**2+(e22/self.delta_d)**2#+(e23/self.delta_m)**2
            de3 = (e31/self.delta_d)**2+(e32/self.delta_d)**2#+(e33/self.delta_m)**2
            de4 = (e41/self.delta_d)**2+(e42/self.delta_d)**2#+(e43/self.delta_m)**2
            de5 = (e51/self.delta_d)**2+(e52/self.delta_d)**2#+(e53/self.delta_m)**2
            de6 = (e61/self.delta_d)**2+(e62/self.delta_d)**2#+(e63/self.delta_m)**2
            de7 = (e71/self.delta_d)**2+(e72/self.delta_d)**2#+(e73/self.delta_m)**2
            de8 = (e81/self.delta_d)**2+(e82/self.delta_d)**2#+(e83/self.delta_m)**2
            de = sqrt(min([de1, de2, de3, de4, de5, de6, de7, de8]))
            gx = min([e11, e21, e31, e41, e51, e61, e71, e81])
            gy = min([e12, e22, e32, e42, e52, e62, e72, e82])
            gz = min([e13, e23, e33, e43, e53, e63, e73, e83])
            dee = (gx/self.delta_d)**2 + (gy/self.delta_d)**2# + (gz/self.delta_m)**2
            #print sqrt(dee), de
            if (x0 + i < 0 or y0 + j < 0 or
                x0 + i >= width or y0 + j >= height or
                z0 + k < zrange[0] or z0 + k >= zrange[1]):
              value2 = MaskCode.Valid
            else:
              if (de <= 1.0):
                value2 = MaskCode.Valid | MaskCode.Foreground
              else:
                value2 = MaskCode.Valid | MaskCode.Background
            new_mask[k,j,i] = value2

      try:
        assert(all(m1 == m2 for m1, m2 in zip(mask, new_mask)))
      except Exception:
        import numpy
        numpy.set_printoptions(threshold=10000)
        diff = (mask == new_mask).as_numpy_array()
        print diff.astype(numpy.int)
        #print mask.as_numpy_array()
        #print new_mask.as_numpy_array()
        #print (new_mask.as_numpy_array()[:,:,:] %2) * (new_mask.as_numpy_array() == 5)
        raise

    # Test passed
    print 'OK'

  def generate_reflections(self, num):
    from random import randint, seed
    from scitbx import matrix
    from dials.array_family import flex
    from dials.algorithms.shoebox import MaskCode
    seed(0)
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

    sigma_b = self.experiment[0].beam.get_sigma_divergence(deg=False)
    sigma_m = self.experiment[0].crystal.get_mosaicity(deg=False)

    rlist = flex.reflection_table()
    rlist['id'] = flex.size_t(len(beam_vector), 0)
    rlist['s1'] = beam_vector
    rlist['panel'] = panel
    rlist['xyzcal.px'] = xyzcal_px
    rlist['xyzcal.mm'] = xyzcal_mm
    rlist['bbox'] = rlist.compute_bbox(self.experiment)
    index = []
    image_size = self.experiment[0].detector[0].get_image_size()
    array_range = self.experiment[0].scan.get_array_range()
    bbox = rlist['bbox']
    for i in range(len(rlist)):
      x0, x1, y0, y1, z0, z1 = bbox[i]
      if (x0 < 0 or x1 > image_size[0] or
          y0 < 0 or y1 > image_size[1] or
          z0 < array_range[0] or z1 > array_range[1]):
        index.append(i)
    rlist.del_selected(flex.size_t(index))
    rlist['shoebox'] = flex.shoebox(
      rlist['panel'], rlist['bbox'])
    rlist['shoebox'].allocate_with_value(MaskCode.Valid)
    return rlist

if __name__ == '__main__':
  from dials.test import cd_auto
  with cd_auto(__file__):
    test = Test()
    test.run()

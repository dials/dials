#!/usr/bin/env python
#
#
#  Copyright (C) (2014) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman.
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import absolute_import, division, print_function

from math import sqrt

from scitbx import matrix

class Test(object):

  def __init__(self, test_nave_model = False):

    # Set up experimental models with regular geometry
    from dxtbx.model import BeamFactory
    from dxtbx.model import GoniometerFactory
    from dxtbx.model import DetectorFactory

    # Beam along the Z axis
    self.beam = BeamFactory.make_beam(unit_s0 = matrix.col((0, 0, 1)),
                                       wavelength = 1.0)

    # Goniometer (used only for index generation) along X axis
    self.goniometer = GoniometerFactory.known_axis(matrix.col((1, 0, 0)))

    # Detector fast, slow along X, -Y; beam in the centre, 200 mm distance
    dir1 = matrix.col((1, 0, 0))
    dir2 = matrix.col((0, -1, 0))
    n = matrix.col((0, 0, 1))
    centre = matrix.col((0, 0, 200))
    npx_fast = npx_slow = 1000
    pix_size = 0.2
    origin = centre - (0.5 * npx_fast * pix_size * dir1 +
                       0.5 * npx_slow * pix_size * dir2)
    self.detector = DetectorFactory.make_detector("PAD",
                        dir1, dir2, origin,
                        (pix_size, pix_size),
                        (npx_fast, npx_slow),
                        (0, 1.e6))

    # Cubic 100 A^3 crystal
    a = matrix.col((100, 0, 0))
    b = matrix.col((0, 100, 0))
    c = matrix.col((0, 0, 100))

    if test_nave_model:
      from dxtbx.model import MosaicCrystalSauter2014
      self.crystal = MosaicCrystalSauter2014(a, b, c, space_group_symbol = "P 1")
      self.crystal.set_half_mosaicity_deg(500)
      self.crystal.set_domain_size_ang(0.2)
    else:
      from dxtbx.model import Crystal
      self.crystal = Crystal(a, b, c, space_group_symbol = "P 1")

    # Collect these models in an Experiment (ignoring the goniometer)
    from dxtbx.model.experiment_list import Experiment
    self.experiment = Experiment(beam=self.beam, detector=self.detector,
      goniometer=None, scan=None, crystal=self.crystal, imageset=None)

    # Generate some reflections
    self.reflections = self.generate_reflections()

    return

  def generate_reflections(self):
    """Use reeke_model to generate indices of reflections near to the Ewald
    sphere that might be observed on a still image. Build a reflection_table
    of these."""
    from cctbx.sgtbx import space_group_info

    space_group_type = space_group_info("P 1").group().type()

    # create a ReekeIndexGenerator
    UB = self.crystal.get_A()
    axis = self.goniometer.get_rotation_axis()
    s0 = self.beam.get_s0()
    dmin = 1.5
    # use the same UB at the beginning and end - the margin parameter ensures
    # we still have indices close to the Ewald sphere generated
    from dials.algorithms.spot_prediction import ReekeIndexGenerator
    r = ReekeIndexGenerator(UB, UB, space_group_type, axis, s0, dmin=1.5, margin=1)

    # generate indices
    hkl = r.to_array()
    nref = len(hkl)

    # create a reflection table
    from dials.array_family import flex
    table = flex.reflection_table()
    table['flags'] = flex.size_t(nref, 0)
    table['id']    = flex.int(nref, 0)
    table['panel'] = flex.size_t(nref, 0)
    table['miller_index'] = flex.miller_index(hkl)
    table['entering']     = flex.bool(nref, True)
    table['s1']           = flex.vec3_double(nref)
    table['xyzcal.mm']    = flex.vec3_double(nref)
    table['xyzcal.px']    = flex.vec3_double(nref)

    return table

  def run(self):

    # cache objects from the model
    UB = matrix.sqr(self.crystal.get_A())
    s0 = matrix.col(self.beam.get_s0())
    es_radius = s0.length()

    # create the predictor and predict for reflection table
    from dials.algorithms.spot_prediction import StillsReflectionPredictor
    predictor = StillsReflectionPredictor(self.experiment)
    predictor.for_reflection_table(self.reflections, UB)

    # for every reflection, reconstruct relp rotated to the Ewald sphere (vector
    # r) and unrotated relp (vector q), calculate the angle between them and
    # compare with delpsical.rad
    from libtbx.test_utils import approx_equal
    for ref in self.reflections:
      r = matrix.col(ref['s1']) - s0
      q = UB * matrix.col(ref['miller_index'])
      tst_radius = (s0 + q).length()
      sgn = -1 if tst_radius > es_radius else 1
      delpsi = sgn*r.accute_angle(q)
      assert approx_equal(delpsi, ref['delpsical.rad'])


  def spherical_relp(self):

    # cache objects from the model
    UB = matrix.sqr(self.crystal.get_A())
    s0 = matrix.col(self.beam.get_s0())
    es_radius = s0.length()

    # create the predictor and predict for reflection table
    from dials.algorithms.spot_prediction import StillsReflectionPredictor
    predictor = StillsReflectionPredictor(self.experiment, spherical_relp=True)
    predictor.for_reflection_table(self.reflections, UB)

    # for every reflection, reconstruct relp centre q, calculate s1 according
    # to the formula in stills_prediction_nave3.pdf and compare
    from libtbx.test_utils import approx_equal
    for ref in self.reflections:
      q = UB * matrix.col(ref['miller_index'])
      radicand = q.length_sq() + 2.0 * q.dot(s0) + s0.length_sq()
      assert radicand > 0.0
      denom = sqrt(radicand)
      s1 = es_radius * (q + s0) / denom
      assert approx_equal(s1, ref['s1'])



if __name__ == '__main__':

  test = Test()
  test.run()
  test.spherical_relp()

  test = Test(test_nave_model=True)
  test.run()

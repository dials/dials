#!/usr/bin/env python
#
# reciprocal_space.py
#
#  Copyright (C) 2014 Diamond Light Source, James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division


class Simulator(object):
  ''' Class to help with simulation from reciprocal space. '''

  def __init__(self, experiment, sigma_b, sigma_m, n_sigma):
    ''' Initialise with models and parameters. '''
    self.experiment = experiment
    self.sigma_b = sigma_b
    self.sigma_m = sigma_m
    self.n_sigma = n_sigma

  def with_given_intensity(self, N, I, B):
    ''' Generate reflections with a given intensity and background. '''
    from dials.array_family import flex
    return self.with_individual_given_intensity(
      N,
      flex.int(N, I),
      flex.int(N, B))

  def with_random_intensity(self, N, Imax, Bmax):
    ''' Generate reflections with a random intensity and background. '''
    from dials.array_family import flex
    return self.with_individual_given_intensity(
      N,
      flex.random_size_t(N, Imax).as_int(),
      flex.random_size_t(N, Bmax).as_int())

  def with_individual_given_intensity(self, N, I, B):
    ''' Generate reflections with given intensity and background. '''
    from dials.algorithms.shoebox import MaskForeground
    from dials.array_family import flex
    from dials.util.command_line import Command, ProgressBar
    from dials.algorithms.reflection_basis import CoordinateSystem
    import random

    # Check the lengths
    assert(N == len(I))
    assert(N == len(B))

    # Generate some predictions
    refl = self.generate_predictions(N)

    # Calculate the signal
    progress = ProgressBar(title='Calculating signal for %d reflections' % len(refl))
    m2 = self.experiment.goniometer.get_rotation_axis()
    s0 = self.experiment.beam.get_s0()
    s1 = refl['s1']
    phi = refl['xyzcal.mm'].parts()[2]
    bbox = refl['bbox']
    shoebox = refl['shoebox']
    for i in range(len(refl)):
      cs = CoordinateSystem(m2, s0, s1[i], phi[i])
      for j in range(I[i]):
        e = (random.gauss(0, self.sigma_b),
             random.gauss(0, self.sigma_b),
             random.gauss(0, self.sigma_m))
        s1_dash, phi_dash = cs.to_beam_vector_and_rotation_angle(e)
        x, y = self.experiment.detector[0].get_ray_intersection(s1_dash)
        x, y = self.experiment.detector[0].millimeter_to_pixel((x, y))
        z = self.experiment.scan.get_array_index_from_angle(phi_dash, deg=False)
        if (x < bbox[i][0] or x > bbox[i][1] or
            y < bbox[i][2] or y > bbox[i][3] or
            z < bbox[i][4] or z > bbox[i][5]):
          continue
        x = int(x - bbox[i][0])
        y = int(y - bbox[i][2])
        z = int(z - bbox[i][4])
        shoebox[i].data[z,y,x] += 1
      progress.update(100.0 * float(i) / len(refl))
    progress.finished('Calculated signal impacts for %d reflections' % len(refl))

    # Calculate the background
    progress = ProgressBar(title='Calculating background for %d reflections' % len(refl))
    for i in range(len(refl)):
      n = len(shoebox[i].data)
      for k in flex.random_size_t(n * B[i], n):
        shoebox[i].data[k] += 1
      progress.update(100.0 * float(i) / len(refl))
    progress.finished('Calculated background for %d reflections' % len(refl))

    # Save the expected intensity and background
    refl['intensity.sim'] = I
    refl['background.sim'] = B

    # Return the reflections
    return refl

  def generate_predictions(self, N):
    ''' Generate some reflections. '''
    from dials.algorithms.shoebox import MaskForeground
    from dials.array_family import flex
    from dials.util.command_line import Command
    from dials.algorithms.reflection_basis import CoordinateSystem
    from dials.algorithms import filtering
    import random

    # Generate a list of reflections
    refl = flex.reflection_table.from_predictions([self.experiment])

    # Filter by zeta
    zeta = 0.05
    Command.start('Filtering by zeta >= %f' % zeta)
    mask = filtering.by_zeta(
      self.experiment.goniometer,
      self.experiment.beam,
      refl['s1'], zeta)
    refl.del_selected(mask != True)
    Command.end('Filtered %d reflections by zeta >= %f' % (len(refl), zeta))

    # Sample if specified
    index = random.sample(range(len(refl)), N)
    refl = refl.select(flex.size_t(index))

    # Compute the bounding box
    refl.compute_bbox(self.experiment, self.n_sigma, self.sigma_b, self.sigma_m)

    # Create a load of shoeboxes
    Command.start('Creating shoeboxes for %d reflections' % len(refl))
    refl['shoebox'] = flex.shoebox(refl['panel'], refl['bbox'])
    refl['shoebox'].allocate_with_value((1 << 0))
    Command.end('Created shoeboxes for %d reflections' % len(refl))

    # Get the function object to mask the foreground
    Command.start('Masking Foreground for %d reflections' % len(refl))
    mask_foreground = MaskForeground(
      self.experiment.beam,
      self.experiment.detector,
      self.experiment.goniometer,
      self.experiment.scan,
      self.n_sigma * self.sigma_b,
      self.n_sigma * self.sigma_m)

    # Mask the foreground
    mask_foreground(
      refl['shoebox'],
      refl['s1'],
      refl['xyzcal.px'].parts()[2])
    Command.end('Masked foreground for %d reflections' % len(refl))

    # Return the reflections
    return refl


if __name__ == '__main__':

  from math import pi
  from dials.model.experiment.experiment_list import ExperimentListFactory
  experiments = ExperimentListFactory.from_json_file(
    "/home/upc86896/Projects/cctbx/sources/dials_regression/centroid_test_data/experiments.json")
  sigma_b = 0.058 * pi / 180
  sigma_m = 0.157 * pi / 180
  n_sigma = 3

  N = 100
  I = 1000
  B = 10

  simulate = Simulator(experiments[0], sigma_b, sigma_m, n_sigma)
  simulate.with_random_intensity(N, I, B)
#  simulate(experiments[0], sigma_b, sigma_m, n_sigma, N, I, B)

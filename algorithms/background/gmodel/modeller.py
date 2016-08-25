#!/usr/bin/env python
#
# algorithm.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division

from dials.phil import parse

phil_scope = parse('''

  min_count = 5
    .type = int(value_min=1)

  nsigma = 6
    .type = float(value_min=0)

  sigma = 0.5
    .type = float(value_min=0)

  kernel_size = 9
    .type = int

  niter = 10
    .type = int

''')

class Modeller(object):

  def __init__(self, beam, detector,
               min_count=5, nsigma=6, sigma=0.5,
               kernel_size=9, niter=10):
    from dials.algorithms.background.gmodel import PixelFilter

    self.beam = beam
    self.detector = detector

    self.min_count = min_count
    self.nsigma = nsigma
    self.sigma = sigma
    self.kernel_size = kernel_size
    self.niter = niter

    width, height = detector[0].get_image_size()

    self._filter = PixelFilter(height, width)

    self.detector_mask = None

  def add_image(self, frame, image, mask, reflections):

    height, width = image.all()

    _,_,_,_,z0,z1 = reflections['bbox'].parts()
    selection = (z0 <= frame) & (z1 > frame)
    subset = reflections.select(selection)
    sbox_mask = subset['shoebox'].apply_background_mask(frame, 1, (height, width))

    if self.detector_mask is None:
      self.detector_mask = mask

    mask = mask & sbox_mask

    self._filter.add(image, mask)

  def compute(self):
    from dials.algorithms.background.gmodel import FillGaps
    from dials.algorithms.image.fill_holes import simple_fill
    result = self._filter.compute(self.min_count, self.nsigma)

    data = result.data()
    mask = result.mask()

    data = simple_fill(data, mask)

    fill = FillGaps(self.beam, self.detector[0])

    mask = mask.as_1d().as_int()
    mask = mask - (~self.detector_mask).as_1d().as_int()
    mask.reshape(data.accessor())

    fill(data, mask, self.sigma, self.kernel_size, self.niter)

    return data


class Creator(object):

  def __init__(self, experiment, params):

    self.modeller = None
    self.background = None
    self.experiment = experiment
    self.params = params

  def initialized(self):
    return self.modeller is not None

  def finalized(self):
    return self.background is not None

  def initialize(self):
    assert not self.initialized()
    self.modeller = Modeller(
      self.experiment.beam,
      self.experiment.detector,
      self.params.min_count,
      self.params.nsigma)

  def next_image(self, frame, image, mask, reflections):
    assert self.initialized()
    self.modeller.add_image(frame, image, mask, reflections)

  def finalize(self):
    assert self.initialized()
    self.background = self.modeller.compute()
    self.modeller = None
    if self.params.debug.output:
      import cPickle as pickle
      from logging import info
      filename = self.params.debug.filename
      info("Writing background model to %s" % filename)
      with open(filename, "w") as outfile:
        pickle.dump(self.background, outfile, protocol=pickle.HIGHEST_PROTOCOL)

  def compute(self, reflections):
    from dials.algorithms.background.gmodel import Fitter
    from dials.array_family import flex
    assert self.finalized()
    fitter = Fitter(self.background)
    scale = fitter(reflections['shoebox'])
    success = scale >= 0
    mean = flex.double([
      flex.mean(sbox.background)
      for sbox in reflections['shoebox']
    ])
    reflections['background.mean'] = mean
    reflections['background.scale'] = scale
    reflections.set_flags(success != True, reflections.flags.dont_integrate)
    return success


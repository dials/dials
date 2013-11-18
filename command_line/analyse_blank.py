#!/usr/bin/env python
#
# analyse_blank.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division
from dials.util.script import ScriptRunner

class Script(ScriptRunner):
  '''A class for running the script.'''

  def __init__(self):
    '''Initialise the script.'''

    # Initialise the base class
    ScriptRunner.__init__(self)

  def main(self, params, options, args):
    '''Execute the script.'''
    from dxtbx.imageset import ImageSetFactory
    from dials.algorithms.image.filter import MeanAndVarianceFilterMasked
    from scitbx.array_family import flex

    # Check if there are enough arguments
    if len(args) < 1:
      self.config().print_help()
      return

    # Get the sweep
    sweep = ImageSetFactory.new(args)[0]
    image = sweep[0].as_double()
    mask  = image >= 0
    mask  = mask.as_1d().as_int()
    mask.reshape(flex.grid(image.all()))

    # Get the kernel size
    kernel_size = params.spotfinder.threshold.kernel_size
    mins = (2*kernel_size[0]+1) * (2*kernel_size[1]+1)

    # Get the filtered image
    filtered = MeanAndVarianceFilterMasked(image, mask, kernel_size, mins)
    mean = filtered.mean()
    var  = filtered.sample_variance()
    mask = filtered.mask()
    mean = mean * mask.as_double()
    var  = var * mask.as_double()
    diff = flex.abs(mean - var)

    # Print some statistics
    print "Mean abs(mean - var)", flex.mean(diff.as_1d())
    print "Median abs(mean - var)", flex.median(diff.as_1d())

    # Get the range in which to print
    mv = flex.mean_and_variance((image * mask.as_double()).as_1d())
    vmax = mv.mean() + 3 * mv.unweighted_sample_standard_deviation()

    # Display images
    from matplotlib import pylab
    pylab.subplot(1, 3, 1)
    pylab.title('Local mean')
    pylab.imshow(mean.as_numpy_array(), interpolation='none',
        vmin=0, vmax=vmax)
    pylab.subplot(1, 3, 2)
    pylab.title('Local variance')
    pylab.imshow(var.as_numpy_array(), interpolation='none',
        vmin=0, vmax=vmax)
    pylab.subplot(1, 3, 3)
    pylab.title('abs(mean - variance)')
    pylab.imshow(diff.as_numpy_array(), interpolation='none',
        vmin=0, vmax=vmax)
    pylab.show()


if __name__ == '__main__':

  script = Script()
  script.run()

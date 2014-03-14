#
# mosflm_like.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: Luiso & James
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division
class MosflmProfileFitting:
  def __init__(self, nblocks):

    self.nblocks = nblocks
  def __call__(self, sweep, crystal, rlist, reference = None):

    from dials.algorithms.integration import flex_2d_layering_n_integrating
    from dials.algorithms.integration.call_mosflm_2d  import mosflm_caller
    flex_2d_layering_n_integrating(rlist)
    xmax, ymax = sweep.get_detector()[0].get_image_size()
    rlist = mosflm_caller(rlist, xmax, ymax, self.nblocks)

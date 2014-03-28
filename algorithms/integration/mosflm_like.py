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
  '''A class to perform 2D mosflm-like profile fitting'''

  def __init__(self, experiment, nblocks):
    self.experiment = experiment
    self.nblocks = nblocks


  #def __call__(self, sweep, crystal, ref_table, reference = None):
  #def __call__(self, sweep, ref_table, reference = None):
  def __call__(self, ref_table, reference = None):
    from dials.algorithms.integration import flex_2d_layering_n_integrating
    from dials.algorithms.integration.call_mosflm_2d  import mosflm_caller
    flex_2d_layering_n_integrating(ref_table)
    xmax, ymax = self.experiment.detector[0].get_image_size()
    ref_table = mosflm_caller(ref_table, xmax, ymax, self.nblocks)

#
# dials.algorithms.background.null_subtractor.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division
from dials.interfaces.background import BackgroundSubtractionInterface


class NullSubtractor(BackgroundSubtractionInterface):
  ''' The NULL background subtractor '''
  def __init__(self, **kwargs):
    ''' Initialise the algorithm. '''
    pass

  def __call__(self, experiment, reflections):
    ''' Do the background subtraction (i.e. set all the background to 0)

    Params:
        experiment The experiment data
        reflections The reflections to process

    Returns:
        The background subtracted reflection list

    '''
    from dials.algorithms.background import set_shoebox_background_value
    set_shoebox_background_value(reflections['shoebox'], 0)
    return reflections

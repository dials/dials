#
# dials.algorithms.background.xds_subtractor.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division
from dials.interfaces.background import BackgroundSubtractionInterface


class XdsSubtractor(BackgroundSubtractionInterface):
  ''' The XDS background subtractor '''

  def __init__(self, **kwargs):
    ''' Initialise the algorithm. '''
    from dials.algorithms.background import XdsSubtractorAlgorithm

    # Create the algorithm
    self._subtractor = XdsSubtractorAlgorithm(
        min_data=kwargs.get("min_data", 10))

  def __call__(self, sweep, crystal, reflections):
    ''' Do the background subtraction as in XDS

    Params:
        sweep The sweep to process
        crystal The crystal to use
        reflections The reflections to process

    Returns:
        The background subtracted reflection list

    '''
    from dials.util.command_line import Command

    # Do the background subtraction
    Command.start('Calculating reflection background')
    mask = self._subtractor(reflections['shoebox'])
    reflections.del_selected(mask != True)
    Command.end('Calculated {0} background values'.format(len(reflections)))

    # Return the reflections
    return reflections

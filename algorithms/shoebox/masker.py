#
# masker.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division

class MaskerBase(object):
  '''A root class to that does overlap masking'''

  def __init__(self, experiment):
    ''' Initialise the overlap masking algorithm

    Params:
        experiment The experiment data
    '''
    from dials.algorithms.shoebox import MaskOverlapping

    # Construct the overlapping reflection mask
    self.mask_overlapping = MaskOverlapping()


  def __call__(self, reflections, adjacency_list=None):
    ''' Mask the given reflections.

    Params:
        reflections The reflection list
        adjacency_list The adjacency_list (optional)

    Returns:
        The masked reflection list

    '''
    from dials.util.command_line import Command

    # Mask the overlaps if an adjacency list is given
    if adjacency_list:
      Command.start('Masking overlapping reflections')
      self.mask_overlapping(
        reflections['shoebox'],
        reflections['xyzcal.px'],
        adjacency_list)
      Command.end('Masked {0} overlapping reflections'.format(
          len(adjacency_list)))

    # Return the reflections
    return reflections

class Masker3DProfile(MaskerBase):
  '''A class to perform 3D profile masking'''

  def __init__(self, experiment, delta_d, delta_m):
    ''' Initialise the masking algorithms

    Params:
        experiment The experiment data
        delta_d The extent of the reflection in reciprocal space
        delta_m The extent of the reflection in reciprocal space

    '''
    super(Masker3DProfile, self).__init__(experiment)

    from dials.algorithms.shoebox import MaskForeground

    # Construct the foreground pixel mask
    self.mask_foreground = MaskForeground(
        experiment.beam, experiment.detector,
        experiment.goniometer, experiment.scan,
        delta_d, delta_m)

  def __call__(self, reflections, adjacency_list=None):
    ''' Mask the given reflections.

    Params:
        reflections The reflection list
        adjacency_list The adjacency_list (optional)

    Returns:
        The masked reflection list

    '''
    reflections = super(Masker3DProfile, self).__call__(reflections, adjacency_list)

    from dials.util.command_line import Command

    if self.mask_foreground:
      # Mask the foreground region
      Command.start('Masking foreground pixels')
      self.mask_foreground(
        reflections['shoebox'],
        reflections['s1'],
        reflections['xyzcal.px'].parts()[2])
      Command.end('Masked foreground pixels for {0} reflections'.format(
        len(reflections)))

    # Return the reflections
    return reflections

class MaskerEmpirical(MaskerBase):
  '''A class to perform empirical masking'''

  def __init__(self, experiment, reference):
    ''' Initialise the masking algorithms

    Params:
        experiment The experiment data

    '''
    super(MaskerEmpirical, self).__init__(experiment)

    from dials.algorithms.shoebox import MaskEmpirical

    # Construct the foreground pixel mask
    self.mask_empirical = MaskEmpirical(reference)

  def __call__(self, reflections, adjacency_list=None):
    ''' Mask the given reflections.

    Params:
        reflections The reflection list
        adjacency_list The adjacency_list (optional)

    Returns:
        The masked reflection list

    '''
    reflections = super(MaskerEmpirical, self).__call__(reflections, adjacency_list)

    from dials.util.command_line import Command

    if self.mask_empirical:
      # Mask the foreground region
      Command.start('Masking foreground pixels')
      self.mask_empirical(reflections)
      Command.end('Masked foreground pixels for {0} reflections'.format(
        len(reflections)))

    # Return the reflections
    return reflections

from __future__ import absolute_import, division

class NullBackgroundExt(object):
  ''' An extension class implementing Null background subtraction. '''

  name = 'null'

  def __init__(self, params, experiments):
    ''' Initialise the algorithm.

    :param params: The input phil parameters
    :param experiments: The experiment list

    '''
    pass

  def compute_background(self, reflections, image_volume=None):
    '''
    Compute the background.

    :param reflections: The list of reflections

    '''
    from dials.algorithms.background import set_shoebox_background_value
    from dials.array_family import flex
    set_shoebox_background_value(reflections['shoebox'], 0)
    return flex.bool(len(reflections), True)

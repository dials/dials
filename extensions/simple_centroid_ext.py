from __future__ import absolute_import, division

class SimpleCentroidExt(object):
  ''' An extension class implementing a simple centroid algorithm. '''

  name = 'simple'

  default = True

  def __init__(self, params, experiments):
    ''' Initialise the algorithm.

    :param params: The input phil parameters
    :param experiments: The experiment list

    '''
    self.experiments = experiments

  def compute_centroid(self, reflections, image_volume=None):
    '''
    Compute the centroid.

    :param reflections: The list of reflections

    '''
    from dials.algorithms.centroid.simple.algorithm import Algorithm
    algorithm = Algorithm(self.experiments)
    return algorithm(reflections, image_volume=image_volume)

#!/usr/bin/env python
#
# standard_interfaces.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from dials.framework import plugin

class SpotFinderThreshold(object):
    '''The threshold for the spotfinding algorithm.'''

    __metaclass__ = plugin.Interface

    name = 'threshold'

    @plugin.abstractmethod
    def __call__(self, image):
        pass


def SpotFinderFilter(object)
    '''The filtering for the spotfinding algorithm.'''

    __metaclass__ = plugin.Interface

    name = 'filter'

    @plugin.abstractmethod
    def __call__(self, flags, **kwargs):
        pass


def RefinementParameterisation(object):
    '''The refinement parameterisation.'''

    name = 'parameterisation'


def IntegrationBoundingBox(object):
    '''The bounding box algorithm.'''

    __metaclass__ = plugin.Interface

    name = 'bounding_box'

    @plugin.abstractmethod
    def __call__(self):
        pass


def IntegrationCentroids(object):
    '''The centroid algorithm.'''

    __metaclass__ = plugin.Interface

    name = 'centroid'

    @plugin.abstractmethod
    def __call__(self)
        pass


def IntegrationBackground(object):
    '''The background modelling algorithm.'''

    __metaclass__ = plugin.Interface

    name = 'background'

    @plugin.abstractmethod
    def __call__(self):
        pass


def IntegrationIntensity(object):
    '''The intensity algorithm.'''

    __metaclass__ = plugin.Interface

    name = 'intensity'

    @plugin.abstractmethod
    def __call__(self):
        pass

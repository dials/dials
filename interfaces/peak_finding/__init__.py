#
# __init__.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import division

class SpotFinderInterface(object):
    '''An interface specification for spot finding classes.'''

    def __init__(self, **kwargs):
        '''Initialise the algorithm with some parameters.

        Params:
            kwargs Key word arguments

        '''
        raise RuntimeError('Please overload!')

    def __call__(self, sweep):
        '''The main function of the spot finder. Select the pixels from
        the sweep and then group the pixels into spots. Return the data
        in the form of a reflection list.

        Params:
            sweep The sweep object

        Returns:
            The reflection list

        '''
        raise RuntimeError('Please overload!')

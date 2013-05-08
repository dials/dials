#!/usr/bin/env python
#
# integrate.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import division
from dials.util.script import ScriptRunner
from dials.interfaces.integration import IntegrationInterface

class IntegratrationAlgorithm(IntegrationInterface):
    '''Select the desired integration algorithm.'''

    def __init__(self, **kwargs):
        '''Initialise the algorithm.'''

        # Get the parameters
        self._params = kwargs['params']

        # Select the algorithm
        self._select_algorithm(**kwargs)

    def _select_algorithm(self, **kwargs):
        '''Select the integration algorithm.'''
        from dials.algorithms.integration.integrate2d import Integrate2d
        from dials.algorithms.integration.integrate3d import Integrate3d

        if self.params.integration.algorithm == '2d':
            self._algorithm = Integrate2d(**kwargs)
        elif self.params.integration.algorithm == '3d':
            self._algorithm = Integrate3d(**kwargs)

    def __call__(self, reflections):
        '''Perform the integration using the selected algorithm.'''
        return self._algorithm(reflections)


class Script(ScriptRunner):
    '''A class for running the script.'''

    def __init__(self):
        '''Initialise the script.'''

        # Initialise the base class
        ScriptRunner.__init__(self)

    def main(self, params, options, args):
        '''Execute the script.'''

        # Initialise the integration algorithm
        self.integrate = IntegratrationAlgorithm(params=params)

        # Perform the integration
        self.integrate(reflections)


if __name__ == '__main__':

    # Execute the script
    script = Script()
    script.run()

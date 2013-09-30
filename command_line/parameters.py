#!/usr/bin/env python
#
# parameters.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division
from dials.util.script import ScriptRunner

class Script(ScriptRunner):
    '''A class for running the script.'''

    def __init__(self):
        '''Initialise the script.'''

        # Initialise the base class
        ScriptRunner.__init__(self, show_config_option=False)

    def main(self, params, options, args):
        '''Execute the script.'''

        # Show the configuration parameters
        self.config().print_phil()


if __name__ == '__main__':

    script = Script()
    script.run()

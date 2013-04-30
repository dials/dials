#!/usr/bin/env python
#
# dials.util.__init__.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division


class HalError(RuntimeError):
    def __init__(self, string=''):

        # Get the username
        try:
            from getpass import getuser
            username = getuser()
        except Exception:
            username = 'Dave'

        # Put in HAL error text.
        text = 'I\'m sorry {0}. I\'m afraid I can\'t do that. {1}'.format(
            username, string)

        # Init base class
        RuntimeError.__init__(self, text)

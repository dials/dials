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

if __name__ == '__main__':

    from dials.util.options import OptionParser
    parser = OptionParser()
    parser.print_system_phil(attributes_level=1)

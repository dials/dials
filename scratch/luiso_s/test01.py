#!/usr/bin/env cctbx.python
# cctbx_introduction_1.py
#

import sys
import exceptions

# tools to use from CBFlib
import pycbf

# and cctbx
from scitbx import matrix

# This was the known UB matrix, derived from XDS processing.

UB = matrix.sqr([1.0, 0.0, 0.0,
                 0.0, 1.0, 0.0,
                 0.0, 0.0, 1.0])
print UB
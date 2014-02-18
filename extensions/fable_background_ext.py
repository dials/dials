#!/usr/bin/env python
#
# fable_background_ext.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import division

from dials.interfaces import BackgroundIface

class FableBackgroundExt(BackgroundIface):

  name = 'fable'

  def compute_background(self, reflections):
    pass

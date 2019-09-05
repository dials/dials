from __future__ import absolute_import, division, print_function

import boost.python

ext = boost.python.import_ext("dials_model_data_ext")
from dials_model_data_ext import *

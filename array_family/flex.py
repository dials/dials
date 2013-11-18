from __future__ import division
from cctbx.array_family.flex import *
from dials.model import data
from dials_array_family_flex_ext import *

# Set the 'real' type to either float or double
if get_real_type() == "float":
  real = float
elif get_real_type() == "double":
  real = double
else:
  raise TypeError('unknown "real" type')

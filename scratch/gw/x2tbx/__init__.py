from __future__ import division
try:
  import boost.python
except Exception:
  ext = None
else:
  ext = boost.python.import_ext("x2tbx_ext", optional = False)

if not ext is None:
  from x2tbx_ext import *


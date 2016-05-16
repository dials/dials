from __future__ import division
try:
  import boost.python
except Exception:
  ext = None
else:
  ext = boost.python.import_ext("recviewer_ext", optional = False)

if ext is not None:
  from recviewer_ext import *

from __future__ import annotations

try:
    import boost_adaptbx.boost.python
except Exception:
    ext = None
else:
    ext = boost_adaptbx.boost.python.import_ext("recviewer_ext", optional=False)

if ext is not None:
    from recviewer_ext import *

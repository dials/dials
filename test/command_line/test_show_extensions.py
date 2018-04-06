from __future__ import absolute_import, division, print_function
from libtbx import easy_run

def test_show_extensions():
    # Call dials.merge_reflection_lists
  easy_run.fully_buffered([
    'dev.dials.show_extensions',
  ]).raise_if_errors()

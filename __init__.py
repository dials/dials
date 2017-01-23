from __future__ import absolute_import, division
import dials.extensions

import logging
try:
  from logging import NullHandler
except ImportError:
  class NullHandler(logging.Handler):
    def emit(self, record):
      pass
logging.getLogger(__name__).addHandler(NullHandler())

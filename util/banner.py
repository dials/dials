from __future__ import absolute_import, division, print_function

import logging
import os

banner = 'DIALS (2018) Acta Cryst. D74, 85-97. https://doi.org/10.1107/S2059798317017235'

def print_banner():
  if os.getenv('DIALS_NOBANNER'):
    return
  d = logging.getLogger('dials')
  d.info(banner)
  logging_to_stdout = any(map(lambda h: isinstance(h, logging.StreamHandler), d.handlers))
  if not logging_to_stdout:
    print(banner)

print_banner()

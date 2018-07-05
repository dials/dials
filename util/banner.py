from __future__ import absolute_import, division, print_function

import atexit
import logging
import os
import sys
import time

super_annoying_banner = '''
############################################################################
#                                                                          #
# The DIALS reference paper has now been published:                        #
#                                                                          #
# DIALS: implementation and evaluation of a new integration package (2018) #
# https://doi.org/10.1107/S2059798317017235        Acta Cryst. D74, 85-97. #
#                                                                          #
############################################################################

'''

less_annoying_banner = 'DIALS (2018) Acta Cryst. D74, 85-97. https://doi.org/10.1107/S2059798317017235'

be_super_annoying = (time.time() < 1533081599 or \
                     '1.14-' in os.getenv('PHENIX_VERSION', '')) # Become less annoying at end of July 2018
if be_super_annoying:
  banner = super_annoying_banner
else:
  banner = less_annoying_banner

def print_banner():
  if os.getenv('DIALS_NOBANNER'):
    return
  d = logging.getLogger('dials')
  d.info(banner)
  logging_to_stdout = any(map(lambda h: isinstance(h, logging.StreamHandler), d.handlers))
  if not logging_to_stdout:
    print(banner)

class ExitHooks(object):
  def __init__(self):
    self.exit_code = None
    self.exception = None

  def hook(self):
    self._orig_exit = sys.exit
    sys.exit = self.exit
    self._orig_excepthook = sys.excepthook
    sys.excepthook = self.exc_handler

  def exit(self, code=0):
    self.exit_code = code
    self._orig_exit(code)

  def exc_handler(self, *args, **kwargs):
    self.exception = True
    return self._orig_excepthook(*args, **kwargs)
hooks = ExitHooks()
hooks.hook()

print_banner()

start = time.time()
@atexit.register
def exit_with_banner():
  if not be_super_annoying: return # no end banner if not annoying
  if hooks.exception is not None: return # no banner on crash
  if hooks.exit_code not in (0, None): return # no banner on ragequit
  if time.time() <= start + 10: return # no banner if we are quick
  print_banner()

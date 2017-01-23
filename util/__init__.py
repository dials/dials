#!/usr/bin/env python
#
# dials.util.__init__.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division
from libtbx.utils import Sorry

from scitbx.array_family import flex # import dependency
from dials_util_ext import *

def debug_console():
  '''Start python console at the current code point.'''
  import sys

  # use exception trick to pick up the current frame
  try:
    raise None
  except: # intentional
    frame = sys.exc_info()[2].tb_frame.f_back

  # evaluate commands in current namespace
  namespace = frame.f_globals.copy()
  namespace.update(frame.f_locals)

  try:
    # Start IPython console if IPython is available.
    from IPython import embed
    embed(user_ns=namespace)
  except ImportError:
    # Otherwise use basic python console
    import code
    code.interact(banner='='*80, local=namespace)

class UsefulError(RuntimeError):
  '''Error message to direct user to report to dials developers.'''

  def __init__(self, message=''):
    if message:
      text = 'Error: "%s" - please report this to dials-support@lists.sourceforge.net' % message
    else:
      text = 'An error has occurred with no message'

    RuntimeError.__init__(self, text)

def usefulraiser(e):
  ''' Function to re-raise an exception with a useful message. '''

  text = 'Please report this error to dials-support@lists.sourceforge.net:'

  if len(e.args) == 0:
    e.args = (text,)
  elif issubclass(e.__class__, Sorry):
    raise
  elif len(e.args) == 1:
    e.args = (text + ' ' + str(e.args[0]),)
  else:
    e.args = (text,) + e.args

  raise

class HalError(RuntimeError):
  def __init__(self, string=''):

    # Get the username
    try:
      from getpass import getuser
      username = getuser()
    except Exception:
      username = 'Dave'

    # Put in HAL error text.
    text = 'I\'m sorry {0}. I\'m afraid I can\'t do that. {1}'.format(
        username, string)

    # Init base class
    RuntimeError.__init__(self, text)

def halraiser(e):
  ''' Function to re-raise an exception with a Hal message. '''

  # Get the username
  try:
    from getpass import getuser
    username = getuser()
  except Exception:
    username = 'Humanoid'

  # Put in HAL error text.
  text = 'I\'m sorry {0}. I\'m afraid I can\'t do that.'.format(username)

  # Append to exception
  if len(e.args) == 0:
    e.args = (text,)
  elif len(e.args) == 1:
    e.args = (text + ' ' + str(e.args[0]),)
  else:
    e.args = (text,) + e.args

  # Reraise the exception
  raise

# clobber existing definitions
HalError = UsefulError
halraiser = usefulraiser

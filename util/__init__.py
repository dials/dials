from __future__ import absolute_import, division, print_function

def debug_console():
  '''Start python console at the current code point.'''
  import sys

  # use exception trick to pick up the current frame
  try:
    raise RuntimeError()
  except RuntimeError:
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

def halraiser(e):
  ''' Function to re-raise an exception with a useful message. '''
  from libtbx.utils import Sorry

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

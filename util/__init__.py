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

# Add the following names to namespace for compatibility reasons.
# Will address relocating them into proper place soonishly. SCRATCH-95
#
def _make_dials_util_ext_redirection(name):
  def dials_util_ext_redirector(*args, **kwargs):
    import dials_util_ext
#   import sys
#   try:
#     raise RuntimeError()
#   except RuntimeError:
#     frame = sys.exc_info()[2].tb_frame.f_back
#   print("DeprecationWarning: {file}:{line} imported method {name} from dials.util rather than from dials_util_ext".format(name=name, file=frame.f_code.co_filename, line=frame.f_lineno))
    return getattr(dials_util_ext, name)(*args, **kwargs)
  return dials_util_ext_redirector
ResolutionMaskGenerator = _make_dials_util_ext_redirection('ResolutionMaskGenerator')
is_inside_polygon = _make_dials_util_ext_redirection('is_inside_polygon')
mask_untrusted_circle = _make_dials_util_ext_redirection('mask_untrusted_circle')
mask_untrusted_polygon = _make_dials_util_ext_redirection('mask_untrusted_polygon')
mask_untrusted_rectangle = _make_dials_util_ext_redirection('mask_untrusted_rectangle')
mask_untrusted_resolution_range = _make_dials_util_ext_redirection('mask_untrusted_resolution_range')
scale_down_array = _make_dials_util_ext_redirection('scale_down_array')

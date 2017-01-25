# LIBTBX_SET_DISPATCHER_NAME dev.dials.remove_absent_reflections

from __future__ import absolute_import, division

from iotbx import mtz

def remove_absent_reflections(hklin, hklout):
  m = mtz.object(hklin)
  s = m.space_group()
  mi = m.extract_miller_indices()

  r = []

  for j, i in enumerate(mi):
    if s.is_sys_absent(i):
      r.append(j)

  for j in reversed(r):
    m.delete_reflection(j)

  m.write(hklout)
  return len(r)

if __name__ == '__main__':
  import os
  import sys
  from libtbx.utils import Sorry

  if len(sys.argv) != 3:
    raise Sorry('%s hklin hklout' % sys.argv[0])

  hklin = sys.argv[1]
  hklout = sys.argv[2]
  if not os.paths.exists(hklin):
    raise Sorry('%s does not exist' % hklin)
  if not os.paths.exists(hklout):
    raise Sorry('%s does not exist' % hklout)

  print 'Removed %d absent reflections' % remove_absent_reflections(
    hklin, hklout)

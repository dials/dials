# LIBTBX_SET_DISPATCHER_NAME dev.dials.pseudo_event_data

from __future__ import absolute_import, division, print_function

def array_to_events(image, n):
  import random
  nslow, nfast = image.focus()

  for s in range(nslow):
    for f in range(nfast):
      if image[(s, f)] <= 0:
        continue
      for j in range(image[(s, f)]):
        print(f, s, n + random.random())

if __name__ == '__main__':
  import sys
  from dxtbx import load
  if len(sys.argv) < 2:
    raise RuntimeError('%s image001.cbf ...' % sys.argv[0])
  for j, arg in enumerate(sys.argv[1:]):
    image = load(arg)
    data = image.get_raw_data()
    array_to_events(data, j)

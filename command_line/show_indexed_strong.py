from __future__ import absolute_import, division

# LIBTBX_SET_DISPATCHER_NAME dev.dials.show_indexed_strong

def show_indexed_strong(indexed_data):

  assert('miller_index' in indexed_data)
  assert('xyzobs.px.value' in indexed_data)

  x_px, y_px, z_px = indexed_data['xyzobs.px.value'].parts()
  mi = indexed_data['miller_index']

  batch = flex.floor(z_px).iround()

  for b in range(min(batch), max(batch) + 1):
    sel = batch == b
    print '%d %d %d' % (b, len(batch.select(sel)), len(mi.select(sel)))

  return

def show_strong(strong_data):

  assert('xyzobs.px.value' in strong_data)

  x_px, y_px, z_px = strong_data['xyzobs.px.value'].parts()

  batch = flex.floor(z_px).iround()

  for b in range(min(batch), max(batch) + 1):
    sel = batch == b
    print '%d %d' % (b, len(batch.select(sel)))

  return

if __name__ == '__main__':
  from dials.array_family import flex # import dependency
  import sys
  if len(sys.argv) != 2:
    raise RuntimeError, '%s indexed.pickle' % sys.argv[0]

  import cPickle as pickle

  data = pickle.load(open(sys.argv[1], 'rb'))

  if 'miller_index' in data:
    show_indexed_strong(data)
  else:
    show_strong(data)

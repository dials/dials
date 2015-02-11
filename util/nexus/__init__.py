from __future__ import division

def open(filename, mode='a'):
  import h5py
  handle = h5py.File(filename, mode)
  if 'entry' in handle:
    entry = handle['entry']
    assert(entry.attrs['NX_class'] == 'NXentry')
  else:
    entry = handle.create_group('entry')
    entry.attrs['NX_class'] = 'NXentry'
  return entry

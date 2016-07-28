from __future__ import division

def enhash(e, h, k, l):
  return e * (2 ** 30) + (h+512) * (2 ** 20) + (k+512) * (2 ** 10) + (l+512)

def dehash(hash):
  e = hash // 2 ** 30
  h = (hash % (2 ** 30)) // (2 ** 20) - 512
  k = (hash % (2 ** 20)) // (2 ** 10) - 512
  l = (hash % (2 ** 10)) - 512
  return e, h, k, l

def add_hash(integrated_data):
  '''Add hash = 2^30 * entering + 2^20 * (h+512) + 2^10 * (k+512) + (l+512)
  as new column to reflection table - should be P1 unique for 360 degree
  scans'''

  from logging import info
  from dials.array_family import flex

  integrated_data = integrated_data.select(integrated_data['id'] >= 0)
  assert max(integrated_data['id']) == 0

  h, k, l = integrated_data['miller_index'].as_vec3_double().parts()
  h = h.iround()
  k = k.iround()
  l = l.iround()
  e = integrated_data['entering'].as_int()

  hash = enhash(e, h, k, l)
  integrated_data['hash'] = hash

  return integrated_data

if __name__ == '__main__':
  import sys
  if len(sys.argv) != 2:
    raise RuntimeError, '%s strong.pickle' % sys.argv[0]

  import cPickle as pickle
  from dials.array_family import flex

  integrated_data = pickle.load(open(sys.argv[1], 'rb'))
  integrated_data = add_hash(integrated_data)
  hash = integrated_data['hash']
  print flex.min(hash), flex.max(hash)

  for h in hash:
    sel = hash == h
    if sel.count(True) > 1:
      print 'Duplicates:', dehash(h)
      isel = sel.iselection()
      for i in isel:
        print integrated_data['miller_index'][i], integrated_data['entering'][i]

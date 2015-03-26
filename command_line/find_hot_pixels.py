from __future__ import division

def print_refl(refl):

  sbox = refl['shoebox'].data
  mask = refl['shoebox'].mask
  ndim = sbox.nd()
  dims = sbox.focus()

  s = refl['shoebox']

  if s.xsize() != 1 or s.ysize() != 1:
    return

  print s.xoffset(), s.yoffset()

  return


if __name__ == '__main__':
  import cPickle as pickle
  import sys
  from dials.array_family import flex # import dependency
  from libtbx.utils import Abort

  if len(sys.argv) != 2:
    raise Abort('exactly 1 reflection table must be specified')

  try:
    table = pickle.load(open(sys.argv[1]))
  except Exception:
    raise Abort('Error loading reflection table')


  for i in range(len(table)):
    print_refl(table[i])

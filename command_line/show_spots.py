# LIBTBX_SET_DISPATCHER_NAME dev.dials.show_spots

from __future__ import division

def show_spots(strong_spots):

  import math

  try:
    x, y, z = zip(*strong_spots['xyzobs.px.value'])
    vx, vy, vz = zip(*strong_spots['xyzobs.px.variance'])
  except RuntimeError as e:
    # convert RuntimeError into more appropriate exception
    raise KeyError(e.message)

  dx = flex.sqrt(flex.double(vx))
  dy = flex.sqrt(flex.double(vy))
  dz = flex.sqrt(flex.double(vz))

  mdx = sum(dx) / len(dx)
  vdx = sum([(v - mdx) ** 2 for v in dx]) / len(dx)

  mdy = sum(dy) / len(dy)
  vdy = sum([(v - mdy) ** 2 for v in dy]) / len(dy)

  mdz = sum(dz) / len(dz)
  vdz = sum([(v - mdz) ** 2 for v in dz]) / len(dz)

  for j in range(len(strong_spots)):
    print '%8.2f %8.2f %8.2f %8.4f %8.4f %8.4f' % (x[j], y[j], z[j],
                                                   dx[j], dy[j], dz[j])

  print '<dX>: %.4f %.4f' % (mdx, math.sqrt(vdx))
  print '<dY>: %.4f %.4f' % (mdy, math.sqrt(vdy))
  print '<dZ>: %.4f %.4f' % (mdz, math.sqrt(vdz))

if __name__ == '__main__':
  import sys
  from libtbx.utils import Sorry
  if len(sys.argv) != 2:
    raise RuntimeError, '%s strong.pickle'

  import cPickle as pickle
  from dials.array_family import flex

  strong_spots = pickle.load(open(sys.argv[1], 'rb'))
  try:
    show_spots(strong_spots)
  except KeyError:
    raise Sorry("{0} does not contain pixel centroid data".format(sys.argv[1]))

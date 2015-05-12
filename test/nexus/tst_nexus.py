
from __future__ import division

if __name__ == '__main__':
  from dials.nexus import load

  filename = '/home/upc86896/Data/NXmx/example.nxs'

  nxmx = load(filename)

  print nxmx[0].title
  print nxmx[0].instrument[0].detector[0].sensor_material
  print len(nxmx[0].instrument[0].detector[0].count_time)
  print nxmx[0].instrument[0].detector[0].module[0].data_size

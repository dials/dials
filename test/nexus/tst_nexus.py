
from __future__ import division

if __name__ == '__main__':
  from dials.nexus import load

  filename = '/home/upc86896/Data/NXmx/example.nxs'

  nxmx = load(filename)

  print nxmx[0].title
  print nxmx[0].instrument

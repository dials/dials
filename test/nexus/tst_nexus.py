
from __future__ import division

if __name__ == '__main__':
  from dials.nexus import Reader

  filename = '/home/upc86896/Data/NXmx/example.nxs'

  reader = Reader(filename)

  nxmx = reader.nxmx()


  print nxmx

  print nxmx.title()
  print nxmx.start_time()
  print nxmx.end_time()
  print nxmx.instrument()
  print nxmx.sample()
  print nxmx.data()

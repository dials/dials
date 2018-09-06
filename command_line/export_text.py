# LIBTBX_SET_DISPATCHER_NAME dev.dials.export_text
from __future__ import absolute_import, division, print_function

from dials.util.export_text import export_text

if __name__ == '__main__':
  import sys
  if len(sys.argv) != 2:
    raise RuntimeError('%s integrated.pickle')

  import six.moves.cPickle as pickle

  with open(sys.argv[1], 'rb') as fh:
    integrated_data = pickle.load(fh)
  export_text(integrated_data)

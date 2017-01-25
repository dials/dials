# LIBTBX_SET_DISPATCHER_NAME dev.dials.export_text
from __future__ import absolute_import, division
from dials.util.export_text import export_text

if __name__ == '__main__':
  import sys
  if len(sys.argv) != 2:
    raise RuntimeError, '%s integrated.pickle'

  import cPickle as pickle

  integrated_data = pickle.load(open(sys.argv[1], 'rb'))
  export_text(integrated_data)

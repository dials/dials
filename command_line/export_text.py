from __future__ import division
from dials.util.export_text import export_text

if __name__ == '__main__':
  import sys
  if len(sys.argv) != 2:
    raise RuntimeError, '%s integrated.pickle'

  import cPickle as pickle
  from dials.array_family import flex

  integrated_data = pickle.load(open(sys.argv[1], 'rb'))
  export_text(integrated_data)

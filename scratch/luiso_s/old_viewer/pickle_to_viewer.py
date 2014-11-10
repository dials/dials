from __future__ import division
from dials.array_family import flex
from dials.viewer.reflection_view import paint_refl

if __name__ == '__main__':
  import cPickle as pickle
  import sys

  table = pickle.load(open(sys.argv[1]))

  ## or
  #table = flex.reflection_table.from_pickle(

  print "num of ref =", len(table)
  tmp_not_needed = '''
  for i in range(len(table)):
    row = table[i]
    print 'Reflection %d' % i
    paint_refl(row, i)
  #'''
  paint_refl(table)

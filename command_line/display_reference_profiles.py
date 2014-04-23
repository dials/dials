#!/usr/bin/env python
#
# display_refernce_profiles.py
#
#  Copyright (C) 2014 Diamond Light Source
#
#  Author: Graeme Winter
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division

def display_reference_profiles(reference_pickle_file):
  '''Display the reference profiles found in the reference_pickle_file generated
  by dials integration using method=fitrs'''

  import cPickle as pickle
  from dials.array_family import flex

  profiles = pickle.load(open(reference_pickle_file))
  central_profile = profiles.profile(5)

  # cast to integer range 0 - 100

  as_integer = 100 * (central_profile / max(central_profile))

  size_z, size_y, size_x = as_integer.focus()

  for k in range(size_z):
    print 'Slice %d' % k
    for j in range(size_y):
      for i in range(size_x):
        print '%4d' % int(as_integer[k, j, i]),
      print ''
    print ''

if __name__ == '__main__':
  import sys
  display_reference_profiles(sys.argv[1])

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

def display_reference_profiles(reference_pickle_file, profile_number,
                               printing=True):
  '''Display the reference profiles found in the reference_pickle_file generated
  by dials integration using method=fitrs'''

  import cPickle as pickle
  from dials.array_family import flex

  profiles = pickle.load(open(reference_pickle_file))
  central_profile = profiles.profile(profile_number)

  # cast to integer range 0 - 100

  as_integer = 100 * (central_profile / max(central_profile))

  size_z, size_y, size_x = as_integer.focus()

  if printing:
    for k in range(size_z):
      print 'Slice %d' % k
      for j in range(size_y):
        for i in range(size_x):
          print '%4d' % int(as_integer[k, j, i]),
        print ''
      print ''

  # now calculate some properties of this profile e.g. the central position and
  # the deviation about this central position

  sum_xi = 0
  sum_yi = 0
  sum_zi = 0

  for k in range(size_z):
    for j in range(size_y):
      for i in range(size_x):
        sum_xi += i * central_profile[k, j, i]
        sum_yi += j * central_profile[k, j, i]
        sum_zi += k * central_profile[k, j, i]

  xi = sum_xi / sum(central_profile)
  yi = sum_yi / sum(central_profile)
  zi = sum_zi / sum(central_profile)

  print 'Moment 1 (zyx): %.2f %.2f %.2f' % (zi, yi, xi)

  sum_xxi = 0
  sum_yyi = 0
  sum_zzi = 0

  for k in range(size_z):
    for j in range(size_y):
      for i in range(size_x):
        sum_xxi += (i - xi) ** 2 * central_profile[k, j, i]
        sum_yyi += (j - yi) ** 2 * central_profile[k, j, i]
        sum_zzi += (k - zi) ** 2 * central_profile[k, j, i]

  import math

  xxi = math.sqrt(sum_xxi / sum(central_profile) - 1)
  yyi = math.sqrt(sum_yyi / sum(central_profile) - 1)
  zzi = math.sqrt(sum_zzi / sum(central_profile) - 1)

  print 'Moment 2 (zyx): %.2f %.2f %.2f' % (zzi, yyi, xxi)

if __name__ == '__main__':
  import sys

  if len(sys.argv) > 2:
    profile_number = int(sys.argv[2])
  else:
    profile_number = 5

  display_reference_profiles(sys.argv[1], profile_number)

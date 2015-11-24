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

# LIBTBX_SET_DISPATCHER_NAME dev.dials.display_reference_profiles

from __future__ import division
import math

def display_reference_profiles(reference_pickle_file, profile_number,
                               printing=True, plot_name=None):
  '''Display the reference profiles found in the reference_pickle_file generated
  by dials integration using intensity.algorithm=fitrs'''

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

  if plot_name is not None:
    import matplotlib
    # http://matplotlib.org/faq/howto_faq.html#generate-images-without-having-a-window-appear
    matplotlib.use('Agg') # use a non-interactive backend

    from matplotlib import pyplot

    # try and make a square-ish grid
    ncols = size_z//int(round(math.sqrt(size_z)))
    nrows = size_z//ncols
    if (nrows*ncols) < size_z:
      nrows += 1

    fig, axes = pyplot.subplots(nrows=nrows, ncols=ncols,
                                sharex=True, sharey=True)

    profile = flex.sqrt(as_integer).as_numpy_array()
    max_value = profile.max()
    for k in range(size_z):
      ax = axes.flat[k]
      slice_z = profile[k,:,:]
      im = ax.imshow(slice_z, interpolation='nearest',
                     vmin=0.0, vmax=max_value,
                     cmap='gray_r')
      ax.set_aspect('equal')

    # hide the unused positions in the subplot grid
    for i, ax in enumerate(axes.flat):
      if i >= size_z:
        ax.axis('off')

    # force the subplots to be square
    pyplot.setp(axes.flat, aspect=1.0, adjustable='box-forced')
    cax, kw = matplotlib.colorbar.make_axes([ax for ax in axes.flat])
    fig.colorbar(im, cax=cax, **kw)
    #pyplot.show()
    pyplot.savefig(plot_name, dpi=300)


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

  _xxi = sum_xxi / sum(central_profile)
  _yyi = sum_yyi / sum(central_profile)
  _zzi = sum_zzi / sum(central_profile)
  if _xxi > 1:
    xxi = math.sqrt(_xxi - 1)
  else:
    xxi = 0.0
  if _yyi > 1:
    yyi = math.sqrt(_yyi - 1)
  else:
    yyi = 0.0
  if _zzi > 1:
    zzi = math.sqrt(_zzi - 1)
  else:
    zzi = 0.0

  print 'Moment 2 (zyx): %.2f %.2f %.2f' % (zzi, yyi, xxi)

if __name__ == '__main__':
  import sys

  if len(sys.argv) < 2:
    print "Usage: dials.display_reference_profiles reflections.pickle"
    exit(0)
  if len(sys.argv) > 2:
    profile_number = int(sys.argv[2])
  else:
    profile_number = 5

  display_reference_profiles(sys.argv[1], profile_number,
                             plot_name="profile.png")

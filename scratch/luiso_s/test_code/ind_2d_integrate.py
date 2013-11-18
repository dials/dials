from __future__ import division
# -*- coding: utf-8 -*-
import numpy

def start(data2d, xcoord, ycoord , x_cm, y_cm, x_sigma, y_sigma):

  data2d = numpy.transpose(data2d)
  single_spot_integrate(data2d, xcoord, ycoord , x_cm, y_cm, x_sigma, y_sigma)

def single_spot_integrate(data2d, xcoord, ycoord , x_cm, y_cm, x_sigma, y_sigma):

  n_x = numpy.size(data2d[:, 0:1])
  n_y = numpy.size(data2d[0:1, :])
  print 'n_x =', n_x
  print 'n_y =', n_y

  data2dtmp = data2d
  data2dsmoth = numpy.zeros(n_x * n_y, dtype = float).reshape(n_x, n_y)
  diffdata2d = numpy.zeros(n_x * n_y, dtype = int).reshape(n_x, n_y)
  diffdata2d_ext = numpy.zeros(n_x * n_y, dtype = int).reshape(n_x, n_y)

  ar_size = 8  #                          This three magical variables MUST be
  ext_area = 4  #                         substituted by hard-coded constants
  threshold_shift = 7.0  #                or must be calculated inside the code

  for pos in range(len(xcoord)):
    x_pos = xcoord[pos]
    y_pos = ycoord[pos]

    from_x = int(x_pos) - ar_size
    x_to = int(x_pos + ar_size + 1)
    from_y = int(y_pos - ar_size)
    y_to = int(y_pos + ar_size + 1)

    for times in range(5):
      for x in range(from_x, x_to, 1):
        for y in range(from_y, y_to, 1):
          pscan = float(numpy.sum(data2dtmp[x - 1:x + 2, y - 1:y + 2]) / 9.0)
          data2dsmoth[x, y] = pscan
      data2dtmp = data2dsmoth

    data2dsmoth[from_x:x_to, from_y:y_to] = data2dsmoth[from_x:x_to, from_y:y_to] + threshold_shift

    for x in range(from_x, x_to, 1):
      for y in range(from_y, y_to, 1):
        if data2d[x, y] > data2dsmoth[x, y]:
          diffdata2d[x, y] = 1

    for x in range(from_x, x_to, 1):
      for y in range(from_y, y_to, 1):
        if diffdata2d[x, y] == 1:
          diffdata2d_ext[x - ext_area:x + ext_area + 1, y - ext_area:y + ext_area + 1] = 1

    x_num_sum = 0.0
    y_num_sum = 0.0
    den_sum = 0.0
    for x in range(from_x, x_to, 1):
      for y in range(from_y, y_to, 1):
        if diffdata2d_ext[x, y] == 1:
          x_num_sum = x_num_sum + float(data2d[x, y]) * float(x)
          y_num_sum = y_num_sum + float(data2d[x, y]) * float(y)
          den_sum = den_sum + float(data2d[x, y])
    if den_sum > 0.0:
      x_cm[pos] = x_num_sum / den_sum
      y_cm[pos] = y_num_sum / den_sum
    else:
      print 'den_sum =', den_sum
      x_cm[pos] = -1
      y_cm[pos] = -1

    x_num_sum = 0.0
    y_num_sum = 0.0
    den_sum = 0.0
    for x in range(from_x, x_to, 1):
      for y in range(from_y, y_to, 1):
        if diffdata2d_ext[x, y] == 1:
          x_num_sum = x_num_sum + float(data2d[x, y]) * (float(x) - float(x_cm[pos])) ** 2.0
          y_num_sum = y_num_sum + float(data2d[x, y]) * (float(y) - float(y_cm[pos])) ** 2.0
          den_sum = den_sum + float(data2d[x, y])
    if den_sum > 0.0:
      x_sigma[pos] = numpy.sqrt(x_num_sum / den_sum)
      y_sigma[pos] = numpy.sqrt(y_num_sum / den_sum)
    else:
      print 'den_sum =', den_sum
      x_sigma[pos] = -1
      y_sigma[pos] = -1

  from matplotlib import pyplot as plt
  paintmask = diffdata2d_ext
  paintmask = paintmask + diffdata2d
  for x in range(0, n_x):
    for y in range(0, n_y):
      if data2dsmoth[x, y] >= 15:
        paintmask[x, y] = paintmask[x, y] + 1

  for pos in range(len(xcoord)):
    nint_xcoord = int(xcoord[pos] + .5)
    nint_ycoord = int(ycoord[pos] + .5)
    paintmask[nint_xcoord, nint_ycoord ] = paintmask[nint_xcoord, nint_ycoord] + 1
  for pos in range(len(xcoord)):
    paintmask[int(x_cm[pos] + 0.5), int(y_cm[pos] + 0.5)] = 5
  print "Plotting paintmask"
  plt.imshow(numpy.transpose(paintmask), interpolation = "nearest", origin = 'lower')
  plt.show()

  print "Plotting diffdata2d"
  plt.imshow(numpy.transpose(diffdata2d), interpolation = "nearest", origin = 'lower')
  plt.show()

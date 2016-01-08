
from __future__ import division

def ice_background(j,i, height, width):
  from numpy.random import poisson
  from math import pi, exp, sqrt, ceil

  A = 100

  x = i
  y = j
  xc = width / 2.0
  yc = height / 2.0
  r = sqrt((xc - x)**2 + (yc - y)**2)
  rc = 20
  rs = 1

  g = A *  exp(-(r-rc)**2 / (2*rs**2))
  return poisson(g)


import cv2

if __name__ == '__main__':
  import numpy

  height = 50
  width = 50
  image = numpy.zeros(dtype="uint8", shape=(height, width))

  X = []
  Y = []
  for j in range(height):
    for i in range(width):
      X.append(i)
      Y.append(j)
      image[j,i] = ice_background(j,i,height,width)

  x0 = 5
  x1 = 15
  y0 = 5
  y1 = 15

  x2 = 30
  x3 = 40
  y2 = 30
  y3 = 40

  mask = numpy.zeros(dtype="uint8", shape=(height, width))
  for j in range(y0,y1):
    for i in range(x0,x1):
      image[j,i] = 255
      mask[j,i] =  1
  for j in range(y2,y3):
    for i in range(x2,x3):
      image[j,i] = 255
      mask[j,i] =  1


  # image = image[0:25,0:25]
  # mask = mask[0:25,0:25]

  print image.dtype

  image = cv2.inpaint(image, mask, 3, cv2.INPAINT_NS)

  from matplotlib import pylab
  pylab.imshow(image, interpolation='none')
  pylab.show()

  Z = image.reshape((width*height,))

  # Z = image.as_1d()
  with open("data.txt", "w") as outfile:
    for x, y, z in zip(X,Y,Z):
      outfile.write("%f %f %f\n" % (x,y,z))

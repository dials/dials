
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

if __name__ == '__main__':

  from dials.array_family import flex
  from math import sqrt

  height = 50
  width = 50
  image = flex.double(flex.grid(height,width))

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

  for j in range(y0,y1):
    for i in range(x0,x1):
      image[j,i] = 0

  data = []
  for j in range(height):
    for i in range(width):
      if (i < x0 or i >= x1) and (j < y0 or j >= y1):
        r = sqrt((i+0.5-width/2.0)**2+(j+0.5-height/2.0)**2)
        data.append((r,image[j,i]))


  data = sorted(data, key=lambda x: x[0])
  radius, data = zip(*data)



  from matplotlib import pylab
  pylab.plot(radius, data)
  pylab.show()

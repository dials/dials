
from __future__ import division

def normal_background(j,i, height, width):
  from numpy.random import poisson
  from math import pi, exp, sqrt, ceil

  A = 100

  x = i
  y = j
  xc = width / 2.0
  yc = height / 2.0
  xs = xc / 4.0
  ys = yc / 4.0

  gc = 1.0 / (2*pi*xs*ys)
  g = A *  exp(-(x-xc)**2 / (2*xs**2))*exp(-(y-yc)**2 / (2*ys**2))
  return poisson(g)

if __name__ == '__main__':

  from dials.array_family import flex

  height = 50
  width = 50
  image = flex.double(flex.grid(height,width))

  X = []
  Y = []
  for j in range(height):
    for i in range(width):
      X.append(i)
      Y.append(j)
      image[j,i] = normal_background(j,i,height,width)

  x0 = 15
  x1 = 25
  y0 = 15
  y1 = 25

  x0 = 20
  x1 = 30
  y0 = 20
  y1 = 30

  for j in range(y0,y1):
    for i in range(x0,x1):
      image[j,i] = 0

  for k in range(100):
    new_image = flex.double(flex.grid(y1-y0,x1-x0))
    for j in range(y0,y1):
      for i in range(x0,x1):
        points = [
          image[j-1,i-1],
          image[j-1,i],
          image[j-1,i+1],
          image[j,i-1],
          image[j,i+1],
          image[j+1,i-1],
          image[j+1,i],
          image[j+1,i+1]]
        j1 = j - y0
        i1 = i - x0
        new_image[j1,i1] = sum(points) / len(points)
    image[y0:y1,x0:y1] = new_image

  from matplotlib import pylab
  #pylab.imshow(image.as_numpy_array(), interpolation='none')
  #pylab.show()

  Z = image.as_1d()
  with open("data.txt", "w") as outfile:
    for x, y, z in zip(X,Y,Z):
      outfile.write("%f %f %f\n" % (x,y,z))


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


def update(image, x0, x1, y0, y1):

  del_I = flex.double(flex.grid(height, width))
  for j in range(height):
    for i in range(width):
      if i > 0 and j > 0 and i < width-1 and j < height-1:
        del_X = (image[j,i+1] - image[j,i-1]) / 2
        del_Y = (image[j+1,i] - image[j-1,i]) / 2
      else:
        del_X = 0
        del_Y = 0
      del_I[j,i] = del_X + del_Y

  K = 1
  C = flex.double(flex.grid(height, width))
  for j in range(height):
    for i in range(width):
      C[j,i] = exp(-(abs(del_I[j,i]) / K)**2)


  delta_I = flex.double(flex.grid(height, width))
  for j in range(height):
    for i in range(width):
      if i > 0 and j > 0 and i < width-1 and j < height-1:
        delta_X = (image[j,i+1] - 2*image[j,i] + image[j,i-1])
        delta_Y = (image[j+1,i] - 2*image[j,i] +  image[j-1,i])
      else:
        delta_X = 0
        delta_Y = 0
      delta_I[j,i] = delta_X + delta_Y

  del_C = flex.double(flex.grid(height, width))
  for j in range(height):
    for i in range(width):
      if i > 0 and j > 0 and i < width-1 and j < height-1:
        del_X = (C[j,i+1] - C[j,i-1]) / 2
        del_Y = (C[j+1,i] - C[j-1,i]) / 2
      else:
        del_X = 0
        del_Y = 0
      del_C[j,i] = del_X + del_Y

  difference = del_C * del_I + C * delta_I
  print flex.max(difference), flex.min(difference)

  image = image + difference
  #image[y0:y1,x0:x1] = image[y0:y1,x0:x1] + difference[y0:y1,x0:x1]

  from matplotlib import pylab
  #pylab.imshow(difference.as_numpy_array())
  #pylab.show()

  return image

if __name__ == '__main__':

  from dials.array_family import flex
  from math import exp

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

  for i in range(10):
    image = update(image, x0, x1, y0, 1)


  # for k in range(100):
  #   new_image = flex.double(flex.grid(y1-y0,x1-x0))
  #   for j in range(y0,y1):
  #     for i in range(x0,x1):
  #       points = [
  #         image[j-1,i-1],
  #         image[j-1,i],
  #         image[j-1,i+1],
  #         image[j,i-1],
  #         image[j,i+1],
  #         image[j+1,i-1],
  #         image[j+1,i],
  #         image[j+1,i+1]]
  #       j1 = j - y0
  #       i1 = i - x0
  #       new_image[j1,i1] = sum(points) / len(points)
  #   image[y0:y1,x0:y1] = new_image

  from matplotlib import pylab
  pylab.imshow(image.as_numpy_array(), interpolation='none')
  pylab.show()

  Z = image.as_1d()
  with open("data.txt", "w") as outfile:
    for x, y, z in zip(X,Y,Z):
      outfile.write("%f %f %f\n" % (x,y,z))

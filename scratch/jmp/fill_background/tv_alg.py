
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

def del_op(A):
  from dials.array_family import flex

  del_A = flex.vec2_double(A.accessor())

  for j in range(A.all()[0]):
    for i in range(A.all()[1]):
      del_x = 0
      del_y = 0
      if i > 0:
        del_x = A[j,i] - A[j,i-1]
      if j < 0:
        del_y = A[j,i] - A[j-1,i]
      del_A[j,i] = (del_y, del_x)

  return del_A

def delta_op(A):
  from dials.array_family import flex

  delta_A = flex.double(A.accessor())

  for j in range(A.all()[0]):
    for i in range(A.all()[1]):
      if j > 0 and j < A.all()[0]-1 and i > 0 and i < A.all()[1]-1:
        delta_A[j,i] = sum([
          -4 * A[j,i],
          A[j,i-1],
          A[j,i+1],
          A[j-1,i],
          A[j+1,i]])

  return delta_A

def div_op(A):

  div_A = flex.double(A.accessor(),0)
  for j in range(A.all()[0]):
    for i in range(A.all()[1]):
      if i > 0 and j > 0:
        div_A[j,i] = A[j,i][0]-A[j-1,i][0]+A[j,i][1]-A[j,i-1][1]

  return div_A


def tv_alg(image, mask, tolerance=1e-3, max_iter=10):
  from dials.array_family import flex
  from scitbx import matrix

  # Set lambda
  L = flex.double(image.accessor(), 0)
  for j in range(image.all()[0]):
    for i in range(image.all()[1]):
      if mask[j,i]:
        L[j,i] = 100
      else:
        L[j,i] = 0

  gamma = 1
  F = image
  UPREV = flex.double(image.accessor(), 0)
  D = flex.vec2_double(image.accessor(), (0,0))
  B = flex.vec2_double(image.accessor(), (0,0))
  for num_iter in range(max_iter):

    U = UPREV

    del_U = del_op(U)

    # Solve D subproblem
    for j in range(image.all()[0]):
      for i in range(image.all()[1]):
        d = matrix.col(del_U[j,i]) + matrix.col(B[j,i])
        dn = d.length()
        if dn == 0:
          D[j,i] = (0,0)
        else:
          D[j,i] = (d / dn) * max([dn-1.0/gamma, 0])

    # Solve U subproblem
    div_db = div_op(D - B)
    delta_U = delta_op(U)

    UCURR = flex.double(image.accessor(), 0)
    for j in range(image.all()[0]):
      for i in range(image.all()[1]):
        if L[j,i] == 0:
          UCURR[j,i] = 0
        else:
          const = gamma / L[j,i]
          UCURR[j,i] = F[j,i] - const*div_db[j,i]+ const*delta_U[j,i]

    # Do the update
    B = B + del_U - D

    # Break if reached convergence
    from math import sqrt
    l2norm = sqrt(sum([u**2 for u in (UCURR - UPREV)]))
    print l2norm
    if l2norm < tolerance:
      break

    UPREV = UCURR
    from matplotlib import pylab
    #pylab.imshow(UCURR.as_numpy_array())
    #pylab.show()

  return UCURR

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
      image[j,i] = ice_background(j,i,height,width)

  x0 = 5
  x1 = 15
  y0 = 5
  y1 = 15


  mask = flex.bool(flex.grid(height, width),True)
  for j in range(y0,y1):
    for i in range(x0,x1):
      image[j,i] = 0
      mask[j,i] =  False

  image = tv_alg(image, mask)

  from matplotlib import pylab
  pylab.imshow(image.as_numpy_array(), interpolation='none')
  pylab.show()

  Z = image.as_1d()
  with open("data.txt", "w") as outfile:
    for x, y, z in zip(X,Y,Z):
      outfile.write("%f %f %f\n" % (x,y,z))

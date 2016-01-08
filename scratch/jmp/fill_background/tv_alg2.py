
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
      if j > 0:
        del_y = A[j,i] - A[j-1,i]
      del_A[j,i] = (del_y, del_x)

  return del_A

def derivative(A, B):
  return A - B

def tv_alg(image, mask, tolerance=1e-3, max_iter=10):
  from scipy.optimize import fmin_bfgs, fmin_l_bfgs_b, fmin_cobyla, fmin_tnc
  from dials.array_family import flex
  from scitbx import matrix
  from math import sqrt

  F = image
  L = flex.double(mask.accessor())
  for i in range(len(mask)):
    if mask[i]:
      L[i] = 5
    else:
      L[i] = 0
  U = flex.double(mask.accessor(), 0)
  for i in range(len(U)):
    U[i] = F[i]

  def func(U):

    U1 = flex.double(F.accessor())
    for i in range(len(F)):
      U1[i] = U[i]
    U = U1

    sum1 = 0
    for f, u, l in zip(F, U, L):
      sum1 += l * (f - u)**2
    sum1 *= 0.5

    del_U = del_op(U)
    sum2 = 0
    for du in del_U:
      sum2 += sqrt(du[0]**2 + du[1]**2)

    Y = sum1 + sum2

    FP = fprime(U)
    print Y
    return Y, FP

  def fprime(U):

    U1 = flex.double(F.accessor())
    for i in range(len(F)):
      U1[i] = U[i]
    U = U1

    def calc_c(U, j, i):
      return 2*U[i,j]**2-2*U[i,j]*(U[i-1,j]+U[i,j-1])+(U[i-1,j]**2+U[i,j-1]**2)

    DU = flex.double(F.accessor(), 0)
    for j in range(DU.all()[0]):
      for i in range(DU.all()[1]):
        if i > 0 and j > 0 and i < DU.all()[1]-1 and j < DU.all()[0]-1:
            C1 = calc_c(U, i,   j) # C_{i,j}
            C2 = calc_c(U, i+1, j) # C_{i+1,j}
            C3 = calc_c(U, i, j+1) # C_{i,j+1}
            DU[j,i] = (2*U[j,i] + U[j,i-1] + U[j-1,i]) / sqrt(C1+0.001) + \
                      (U[j,i] - U[j,i+1]) / sqrt(C2+0.001) + \
                      (U[j,i] - U[j+1,i]) / sqrt(C3+0.001) - \
                      L[j,i] * (F[j,i] - U[j,i])
    # print sum(DU)
    return list(DU)

  x0 = U

  min_F = min(F)
  max_F = max(F)
  bounds = [(min_F, max_F) for x in x0]

  #U = fmin_bfgs(func, x0, epsilon=1e-5*sum(image), maxiter=max_iter)
  U, nfeval, rc = fmin_tnc(func, list(x0), bounds=bounds,
               disp=3)


  U1 = flex.double(F.accessor())
  for i in range(len(F)):
    U1[i] = U[i]
  U = U1

  return U





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

  image = image[0:25,0:25]
  mask = mask[0:25,0:25]

  image = tv_alg(image, mask)

  from matplotlib import pylab
  pylab.imshow(image.as_numpy_array(), interpolation='none')
  pylab.show()

  Z = image.as_1d()
  with open("data.txt", "w") as outfile:
    for x, y, z in zip(X,Y,Z):
      outfile.write("%f %f %f\n" % (x,y,z))


from __future__ import division

def glm33(X, Y, B, P):
  from math import exp , sqrt, log, factorial, lgamma
  from scitbx import matrix
  from scipy.stats import poisson

  X = matrix.rec(list(X), X.all())
  Y = matrix.col(list(Y))
  B = matrix.col(list(B))
  P = matrix.col(list(P))

  while True:
    n = X * B
    mu = matrix.col([exp(nn) for nn in n])

    z = [(yy - mui) / mui for yy, mui in zip(Y, mu)]
    w = [pp * mui for pp, mui in zip(P, mu)]

    W = matrix.diag(w)
    Z = matrix.col(z)
    delta = (X.transpose() * W * X).inverse() * X.transpose()*W*z

    relE = sqrt(sum([d*d for d in delta])/max(1e-10, sum([d*d for d in B])))
    B = B + delta

    # print relE
    if relE < 1e-3:
      break

  return B


if __name__ == '__main__':

  from dials.array_family import flex
  from random import uniform
  from numpy.random import poisson, seed
  from math import exp, log
  from scitbx.glmtbx import glm
  seed(0)
  a = 1
  c = 10
  nobs = 100

  X1 = [1] * nobs
  X2 = [i for i in range(nobs)]
  # X1 = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
  # X2 = [0, 0, 1, 2, 4, 5, 4, 3, 2, 0, 0]
  X = flex.double(item for a, b in zip(X1, X2) for item in [a, b])
  X.reshape(flex.grid(nobs, 2))
  # Y = flex.double([poisson(a*X1[i]+c*X2[i],1)[0] for i in range(nobs)])
  Y = flex.double([poisson(c*X1[i]+a*X2[i],1)[0] for i in range(nobs)])
  B = flex.double([0, 0])
  P = flex.double([1] * nobs)
  B1 = glm33(X,Y,B,P)
  print list(B1)
  result = glm(X,Y,B,P, max_iter=1000)
  print list(result.parameters())

  exit(0)
  # for k in range(10):
  #   X = flex.double([item for i in range(100) for item in [1,i-50]])
  #   X.reshape(flex.grid(100, 2))
  #   Y = flex.double([poisson(a*i+c,1)[0] for i in range(100)])
  #   B = flex.double([0, 0])
  #   P = flex.double([1] * len(Y))

  #   B = glm33(X,Y,B,P)
  #   print exp(B[0]), exp(B[1])

  Y2 = []
  Y3 = []
  for i in range(11):
    Y2.append(exp(X[i,0]*B[0]+X[i,1]*B[1]))
    Y3.append(X[i,0]*a+X[i,1]*c)

  from matplotlib import pylab
  pylab.plot(range(11),Y)
  pylab.plot(range(11),Y2)
  pylab.plot(range(11),Y3)
  pylab.show()

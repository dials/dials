
from __future__ import division
def median(x):
  import numpy
  return numpy.median(x)

if __name__ == '__main__':

  from scitbx.glmtbx import robust_glm, glm
  from math import exp, ceil, log
  from dials.array_family import flex
  from numpy.random import poisson
  from time import time

  known = [1, 1]

  X = []
  for i in range(1000):
    X.append([1, i])

  tt = 0

  for it in range(1000):

    # Y = [exp(known[0]*x[0] + known[1]*x[1]) for x in X]

    # Y1 = [int(poisson(exp(known[0]*x[0] + known[1]*x[1]),1)[0]) for x in X]
    Y1 = [int(poisson(known[0]*x[0] + known[1]*x[1],1)[0]) for x in X]
    # Y1 = [int(exp(known[0]*x[0] + known[1]*x[1])) for x in X]

    initial = median(Y1)
    if initial > 0:
      initial = log(initial)

    initial = [initial, 0]

    st = time()
    result = robust_glm(
      flex.double(X),
      flex.double(Y1),
      flex.double(initial),
      max_iter=1000)
    assert result.converged()
    tt += time() - st

    beta = list(result.parameters())

    print beta

    Y2 = [exp(beta[0]*x[0] + beta[1]*x[1]) for x in X]

    P = [1.0 for x in X]

    result = glm(
      flex.double(X),
      flex.double(Y1),
      flex.double(initial),
      flex.double(P),
      max_iter=1000)
    assert result.converged()

    beta = list(result.parameters())

    print beta

    Y3 = [exp(beta[0]*x[0] + beta[1]*x[1]) for x in X]

    m1 = max(Y1)
    m2 = max(Y2)
    m3 = max(Y3)
    m = ceil(max([m1,m2,m3]))

    # from matplotlib import pylab
    # pylab.ylim((0, m))
    # # pylab.plot(Y, label="Y")
    # pylab.plot(Y1, label="Y1")
    # pylab.plot(Y2, label="Y2")
    # pylab.plot(Y3, label="Y3")
    # pylab.legend()
    # pylab.show()

  print "Time: ", time() - st

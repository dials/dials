


from __future__ import division

if __name__ == '__main__':

  from dials.array_family import flex
  # from numpy.random import poisson, seed
  from scitbx.glmtbx import robust_glm
  from math import log
  # seed(0)

  for k in range(10):
    print k
    # a = list(poisson(2, 100))
    a = [3,
 2,
 5,
 1,
 0,
 0,
 7,
 1,
 3,
 3,
 5,
 2,
 2,
 0,
 1,
 3,
 1,
 2,
 1,
 0,
 1,
 1,
 1,
 1,
 5,
 1,
 1,
 1,
 0,
 2,
 1,
 3,
 0,
 2,
 1,
 1,
 2,
 3,
 1,
 4,
 5,
 4,
 1,
 2,
 2,
 2,
 3,
 3,
 2,
 3,
 2,
 0,
 4,
 6,
 2,
 1,
 1,
 2,
 2,
 1,
 2,
 4,
 1,
 2,
 1,
 1,
 2,
 1,
 0,
 0,
 5,
 2,
 2,
 4,
 3,
 3,
 1,
 0,
 1,
 2,
 2,
 1,
 3,
 1,
 5,
 3,
 1,
 0,
 4,
 0,
 4,
 4,
 3,
 0,
 1,
 4,
 1,
 1,
 1,
 6]
    # a[4] = 10
    # a[50] = 100

    X = flex.double([1] * len(a))
    X.reshape(flex.grid(len(a), 1))

    Y = flex.double(a)
    B = flex.double([0])

    result = robust_glm(X, Y, B, family="poisson")
    print list(result.parameters())


  from matplotlib import pylab
  print "MOM1: ", sum(means1) / len(means1)
  print "MOM2: ", sum(means2) / len(means2)
  pylab.plot(means1, color='black')
  pylab.plot(means2, color='blue')
  pylab.show()

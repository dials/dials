

from __future__ import division

def bootstrap(x):
  from random import sample

  n = len(x)

  r = []
  for i in range(len(x)):
    y = x[:i] + x[i+1:]
    r.append(sum(y) / len(y))

  from matplotlib import pylab
  pylab.hist(r, bins=50)
  pylab.show()

  return sum(r) / len(r)




if __name__ == '__main__':

  from dials.array_family import flex
  from random import uniform
  from numpy.random import poisson, seed
  seed(0)
  means1 = []
  means2 = []
  for k in range(100):
    a = list(poisson(0.1, 100))
    a[4] = 100
    means1.append(sum(a)/len(a))
    means2.append(bootstrap(a))

  from matplotlib import pylab
  print "MOM1: ", sum(means1) / len(means1)
  print "MOM2: ", sum(means2) / len(means2)
  pylab.plot(means1, color='black')
  pylab.plot(means2, color='black')
  pylab.show()

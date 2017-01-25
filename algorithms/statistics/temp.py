from __future__ import absolute_import, division

#from scipy.stats.distributions import ksone

#v1 = []
#v2 = []

#x = 0.5
#n = 100
#m = int(floor(n * (1 - x)))
#for j in range(0, m+1):
  #a = x + float(j) / n
  #if a > 0 and a < 1:
    #b = (j - 1) / log(a)
    #c = (n - j) / log(1 - a)
  #elif a <= 0:
    #b = 0
    #c = 0
  #elif a >= 1:
    #b = -100
    #c = 0

  #v1.append(b)
  #v2.append(c)
#from matplotlib import pylab
#pylab.plot(v1)
#pylab.plot(v2)
#pylab.show()

#exit(0)

#def smirnov2(n, e):

  #from math import floor

  #assert(n > 0 and e >= 0.0 and e <= 1.0)

  #if e == 0.0:
    #return 1.0
  #nn = int(floor(n * (1.0 - e)))
  #p = 0.0
  #c = 1.0
  #for v in range(0, nn+1):
    #evn = e + float(v) / n
    #aa = pow(evn, (v - 1))
    #bb = pow(1.0 - evn, n - v)
    #p += c * aa * bb
    #print aa, bb, c, p
    #c *= float(n - v) / (v + 1)
    ##print v, c
  #return p * e

def compute_lz(z):
  from math import sqrt, pi, exp
  s = sum(exp(-(2.0*k-1.0)**2 * pi**2 / (8.0*z**2)) for k in range(1,10000))
  return s * sqrt(2.0*pi) / z

from dials.algorithms.statistics import *
from math import sqrt

xx = []
scdf = []
kcdf = []
x = 0.01

#for n in range(1, 100):
  #s = smirnov_cdf(n, x)
  #k = compute_lz(x) / sqrt(n)
  ##k = kolmogorov_cdf(n, x)
  #xx.append(n)
  #scdf.append(s)
  #kcdf.append(k)

n = 100
for i in range(1, 200):
  x = float(i) / 100.0
  s = 0
#  s = smirnov_cdf(n, x)
  k = 0
  k = compute_lz(x) / sqrt(n)
  xx.append(x)
  scdf.append(s)
  kcdf.append(k)


from matplotlib import pylab
p1 = pylab.plot(xx, scdf, color='r')
p2 = pylab.plot(xx, kcdf, color='b')
pylab.show()

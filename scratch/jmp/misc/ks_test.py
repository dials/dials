

from scipy.special import gammaln

def smirnov(n, e):
  from math import floor, log, exp
  MAXLOG = 2**1022
  assert(n > 0 and e >= 0.0 and e <= 1.0)

  if e == 0.0:
    return 1.0

  nn = int(floor(n * (1.0 - e)))
  p = 0.0
  if n < 1013:
    c = 1.0
    for v in range(0, nn+1):
      evn = e + float(v) / n
      p += c * pow(evn, (v - 1)) * pow(1.0 - evn, n - v)
      c *= float(n - v) / (v + 1)
  else:
    lgamnp1 = gammaln(n + 1)
    for v in range(0, nn+1):
      evn = e + float(v) / n
      omevn = 1.0 - evn
      if abs(omevn) > 0.0:
        t = lgamnp1 - gammaln(v + 1) - gammaln(n - v + 1) + \
          (v - 1) * log(evn) + (n - v) * log(omevn)
        if t > -MAXLOG:
          p += exp(t)
  return p * e


def smirnov2(n, e):

  from math import floor

  assert(n > 0 and e >= 0.0 and e <= 1.0)

  if e == 0.0:
    return 1.0
  nn = int(floor(n * (1.0 - e)))
  p = 0.0
  c = 1.0
  for v in range(0, nn+1):
    evn = e + float(v) / n
    aa = pow(evn, (v - 1))
    bb = pow(1.0 - evn, n - v)
    p += c * aa * bb
    c *= float(n - v) / (v + 1)
    print v, c
  return p * e

def smirnov3(n, e):
  from math import floor, log, exp
  MAXLOG = 2**1022
  nn = int(floor(n * (1.0 - e)))
  p = 0.0
  lgamnp1 = gammaln(n + 1)
  for v in range(0, nn+1):
    evn = e + float(v) / n
    omevn = 1.0 - evn
    if abs(omevn) > 0.0:
      t = lgamnp1 - gammaln(v + 1) - gammaln(n - v + 1) + \
        (v - 1) * log(evn) + (n - v) * log(omevn)
      if t > -MAXLOG:
        p += exp(t)
  return p * e

from scipy.stats.distributions import ksone

Dplus = 0.0002
N = 100

print ksone.sf(Dplus, N)
print smirnov(N, Dplus)

print smirnov2(1300, 0.01)

#for n in range(1000, 1500):
  #print n, smirnov2(n, 0.01)


#N = 10
#v = []
#y = []
#for x in range(0, 0.5):
  #y.append(x / 1000.0)
  #v.append(1.0 - smirnov2(N, x / 1000.0))


#from matplotlib import pylab
#pylab.plot(y, v)
#pylab.show()



#from time import time
#st = time()
#r = smirnov(

#from time import time
#st = time()
#for i in range(1000):
  #smirnov2(1000, Dplus)
#print time() - st

#st = time()
#for i in range(1000):
  #smirnov3(1000, Dplus)
#print time() - st
#for i in range(1, 100):
  #print smirnov2(i, Dplus), smirnov3(i, Dplus), ksone.sf(Dplus, i)

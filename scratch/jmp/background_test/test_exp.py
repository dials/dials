
from __future__ import division

def pmf(mu, k):
  from scipy.stats import poisson
  if k < 0:
    return 0
  return poisson.pmf(k, mu)

def cdf(mu, k):
  from scipy.stats import poisson
  if k < 0:
    return 0
  return poisson.cdf(k, mu)


def expectation(mu, c):
  from math import sqrt, floor
  svar = sqrt(mu)
  j1 = int(floor(mu - c*svar))
  j2 = int(floor(mu + c*svar))
  p1 = pmf(mu, j1)
  p2 = pmf(mu, j2)
  p3 = cdf(mu, j1)
  p4 = pmf(mu, j2+1)
  p5 = cdf(mu, j2+1)
  p6 = 1.0 - p5 + p4
  p7 = pmf(mu, j1-1)
  p8 = pmf(mu, j2-1)
  p9 = cdf(mu, j2-1)
  p10 = p9 - p3 + p1
  epsi1 = c*(p6-p3) + (mu/svar)*(p1-p2)
  epsi2 = c*(p1+p2) + (mu*mu/(svar*svar*svar))*(p10/mu+p7-p1-p8+p2)
  return epsi1, epsi2

if __name__ == '__main__':

  c = 1.345
  e1 = []
  e2 = []
  m = []
  for mu in range(1, 10000):
    m.append(mu / 1000.0)
    e = expectation(mu / 1000.0, c)
    e1.append(e[0])
    e2.append(e[1])
  print e1[:10]
  from matplotlib import pylab
  pylab.plot(m, e1)
  pylab.show()

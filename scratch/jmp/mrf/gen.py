
def mean_and_variance(a):
  sum1 = sum(a)
  sum2 = sum(aa*aa for aa in a)
  mean = sum1 / len(a)
  variance = (sum2 - sum1*sum1/len(a)) / len(a)
  return mean, variance


def m_estimate(x):
  from math import sqrt, erf, exp
  from scipy.stats import poisson
  m0 = sum(x) / len(x)

  while True:
    w = [ 1.0/sqrt(1+0.5*(xx-m0)**2) for xx in x ]
    # w = [ poisson.pmf(xx, m0) for xx in x ]
    m = sum(ww * xx for ww, xx in zip(w, x)) / sum(w)
    if abs(m - m0) < 1e-3:
      break
    m0 = m
  r = [xx - m for xx in x]
  return m, w, r

def glm(x):
  m0 = sum(x) / len(x) + 1
  while True:
    m = 2*abs(m0) - len(x)*m0*m0/sum(x)
    if abs(m - m0) < 1e-3:
      break
    m0 = m
  return m

def glm2(x):
  m0 = sum(x) / len(x) + 1
  while True:
    w = [ 1.0/sqrt(1+0.5*(xx-m0)**2) for xx in x ]
    m = 2*abs(m0) - len(x)*m0*m0/sum(x)
    if abs(m - m0) < 1e-3:
      break
    m0 = m
  return m

def sign(x):
  if x == 0:
    return 0
  elif x < 0:
    return -1
  return 1

def glm33(x, p):
  from math import exp , sqrt, log, factorial, lgamma
  from scitbx import matrix
  from scipy.stats import poisson
  beta0 = 0

  while True:
    n = beta0
    mu = exp(n)


    z = [n + (xx - mu) / mu for xx in x]
    w = [pp * mu for pp in p]

    W = matrix.diag(w)
    Z = matrix.col(z)
    X = matrix.rec([1]*len(z), (len(z),1))
    H = X*(X.transpose() *W* X).inverse() * X.transpose()*W
    H = [H[i+i*len(z)] for i in range(len(z))]
    MSE = sum((xx-mu)**2 for xx in x) / len(x)
    r = [xx - mu for xx in x]
    D = [rr*rr / (1 * MSE) * (hh / (1-hh)**2) for rr, hh in zip(r, H)]
    N = sum(1 for d in D if d > 4 / len(x))
    from scipy.special import iv
    I0 = iv(0,2*mu)
    S = [log(I0) + lgamma(xx+1) - mu - xx*log(mu) for xx in x]
    print "S max: ", max(S)

    # print "DMAX: ", max(D), 4 / len(x), N
    # from matplotlib import pylab
    # pylab.hist(D)
    # pylab.show()
    beta = (X.transpose() * W * X).inverse() * X.transpose()*W*z
    if abs(beta[0] - beta0) < 1e-3:
      break
    beta0 = beta[0]

  return exp(beta[0])

def simple(x):
  from numpy import median
  m = median(x)
  mad = median([abs(xx - m) for xx in x])
  print m, mad
  return m

def glm3(x):
  from scipy.stats import poisson
  p = [1 for xx in x]
  return glm33(x,p)


def mtl(x):
  from math import log, factorial
  from matplotlib import pylab
  n = 5#int(len(x) / 20)
  y = x[:]

  mu0 = sum(y) / len(y)

  while True:

    l = [xx*log(mu0)-mu0-log(factorial(xx)) for xx in x]
    ind = sorted(range(len(l)), key=lambda i: l[i])
    ind = ind[n:]
    y = [x[i] for i in ind]
    mu = sum(y) / len(y)
    if abs(mu - mu0) < 1e-3:
      break
    mu0 = mu
  return mu


if __name__ == '__main__':

  from dials.array_family import flex
  from random import uniform
  from numpy.random import poisson, seed
  seed(0)
  means1 = []
  means3 = []
  means4 = []
  for k in range(100):
    # print "---"
    a = list(poisson(1, 100))

    a[4] = 1000
    # a[5] = 100

    mean_m, weight_m, res_m = m_estimate(a)
    # print mean_m
    # from matplotlib import pylab
    # pylab.plot(res_m)
    # pylab.show()
    means1.append(sum(a)/len(a))
    means3.append(mtl(a))
    mean_m = glm3(a)
    print k, mean_m
    means4.append(mean_m)

  from matplotlib import pylab
  print "MOM1: ", sum(means1) / len(means1)
  print "MOM3: ", sum(means3) / len(means3)
  print "MOM4: ", sum(means4) / len(means4)
  pylab.plot(means1, color='black')
  pylab.plot(means4, color='blue')
  pylab.show()

  # sum1 = sum(a)
  # p = len(a)
  # c_max = 2*(p-1)
  # b_min = (p-2+1/p)
  # print c_max
  # print b_min
  # b = b_min
  # c = c_max

  # l = [(1 - (c / (b + sum1))) * aa for aa in a]

  # pylab.plot(l)
  # pylab.plot(a)
  # pylab.show()

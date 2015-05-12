

from __future__ import division

# def huber(r, c):
#   if (abs(r) < c):
#     return 1
#   else:
#     return c/abs(r)

def sign(x):
  if x == 0:
    return 0
  elif x < 0:
    return -1
  else:
    return 1

def huber(r, c):
  if (abs(r) < c):
    return r
  else:
    return c * sign(r)


def glm(x, p = None):
  from math import exp , sqrt, log, factorial, lgamma
  from scitbx import matrix
  from scipy.stats import poisson
  beta0 = 0

  if p is None:
    p = [1.0] * len(x)

  X = matrix.rec([1]*len(x), (len(x),1))
  c = 1.345

  while True:
    n = beta0
    mu = exp(n)
    z = [n + (xx - mu) / mu for xx in x]
    w = [pp * mu for pp in p]
    r = [(xx - mu) / sqrt(mu) for xx in x]
    w2 = [huber(rr, c) for rr in r]

    W = matrix.diag(w)
    # W2 = matrix.diag(w2)
    Z = matrix.col(z)
    beta = (X.transpose() * W * X).inverse() * X.transpose()*W*z
    print beta
    if abs(beta[0] - beta0) < 1e-3:
      break
    beta0 = beta[0]

  # W12 = matrix.diag([sqrt(ww) for ww in w])
  # H = W12 * X * (X.transpose() * W * X).inverse() * X.transpose() * W12

  mu = exp(n)
  return mu


def glm2(y):
  from math import sqrt, exp, floor, log
  from scipy.stats import poisson
  from scitbx import matrix
  from dials.array_family import flex

  y = flex.double(y)

  x = flex.double([1.0 for yy in y])
  w = flex.double([1.0 for yy in y])

  X = matrix.rec(x, (len(x),1))

  c = 1.345

  beta = matrix.col([0])

  maxiter = 10
  accuracy = 1e-3

  for iter in range(maxiter):

    ni = flex.double([1.0 for xx in x])
    sni = flex.sqrt(ni)

    eta = flex.double(X*beta)

    mu = flex.exp(eta)
    dmu_deta = flex.exp(eta)

    Vmu = mu
    sVF = flex.sqrt(Vmu)
    residP = (y - mu) * sni / sVF

    phi = 1
    sV = sVF * sqrt(phi)
    residPS = residP / sqrt(phi)

    H = flex.floor(mu*ni - c*sni*sV)
    K = flex.floor(mu*ni + c*sni*sV)

    dpH  = flex.double([poisson(mui).pmf(Hi)   for mui,Hi in zip(mu,H)])
    dpH1 = flex.double([poisson(mui).pmf(Hi-1) for mui,Hi in zip(mu,H)])
    dpK  = flex.double([poisson(mui).pmf(Ki)   for mui,Ki in zip(mu,K)])
    dpK1 = flex.double([poisson(mui).pmf(Ki-1) for mui,Ki in zip(mu,K)])
    pHm1 = flex.double([poisson(mui).cdf(Hi-1) for mui,Hi in zip(mu,H)])
    pKm1 = flex.double([poisson(mui).cdf(Ki-1) for mui,Ki in zip(mu,K)])
    pH = pHm1 + dpH # = ppois(H,*)
    pK = pKm1 + dpK # = ppois(K,*)
    E2f = mu*(dpH1 - dpH - dpK1 + dpK) + pKm1 - pHm1
    Epsi = c * (1.0 - pK - pH) + (mu / sV) * (dpH - dpK)
    Epsi2 = c*c * (pH + 1.0 - pK) + E2f
    EpsiS = c*(dpH + dpK) + E2f / sV

    psi = flex.double([huber(rr, c) for rr in residPS])
    cpsi = psi - Epsi
    temp = cpsi * w * sni/sV * dmu_deta
    EEqMat = [0] * len(X)
    for j in range(X.n_rows()):
      for i in range(X.n_columns()):
        k = i + j * X.n_columns()
        EEqMat[k] = X[k] * temp[j]
    EEqMat = matrix.rec(EEqMat, (X.n_rows(), X.n_columns()))
    EEq = []
    for i in range(EEqMat.n_columns()):
      col = []
      for j in range(EEqMat.n_rows()):
        k = i + j * EEqMat.n_columns()
        col.append(EEqMat[k])
      EEq.append(sum(col)/len(col))
    EEq = matrix.col(EEq)
    DiagB = EpsiS /(sni*sV) * w * (ni*dmu_deta)**2
    B = matrix.diag(DiagB)
    H = (X.transpose() * B * X) / len(y)
    dbeta = H.inverse() * EEq
    beta_new = beta + dbeta

    relE = sqrt(sum([d*d for d in dbeta])/max(1e-20, sum([d*d for d in beta])))
    beta = beta_new

    print relE
    if relE < accuracy:
      break


  weights = [min(1, c / abs(r)) for r in residPS]

  eta = flex.double(X*beta)

  mu = flex.exp(eta)
  return mu[0]


if __name__ == '__main__':

  from dials.array_family import flex
  from random import uniform
  from numpy.random import poisson, seed
  seed(0)
  means1 = []
  means2 = []

  for k in range(10):
    print k
    a = list(poisson(0.1, 100))
    print a
    # a[4] = 10
    # a[50] = 100
    means1.append(sum(a)/len(a))
    means2.append(glm2(a))

  from matplotlib import pylab
  print "MOM1: ", sum(means1) / len(means1)
  print "MOM2: ", sum(means2) / len(means2)
  pylab.plot(means1, color='black')
  pylab.plot(means2, color='blue')
  pylab.show()

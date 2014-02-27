



def standard_normal_cdf(x):
  from math import erf, sqrt
  return (1.0 / 2.0) * (erf(x / sqrt(2.0)) + 1.0)


def distance(t, f):
  f2 = [standard_normal_cdf(tt) for tt in t]
  return [ff2 - ff for ff2, ff in zip(f2, f)]

def supremum(x):
  return max([abs(xx) for xx in x])

def kolmogorov_smirnov_test(data):

  data = standardize(data)

  t, f = empirical_distribution_function(data)
  return supremum(distance(t, f))


def empirical_distribution_function(data):
  sorted_data = sorted(data)
  t = sorted_data
  n = float(len(data))
  f = [x / n for x in range(1, len(data)+1)]
  return t, f

def standardize(data):
  from math import sqrt
  sum_x = sum(data)
  sum_x2 = sum([d*d for d in data])
  mean = sum_x / len(data)
  n = len(data)
  sdev = sqrt((sum_x2 - sum_x*sum_x / n) / (n-1))
  return [(d - mean) / sdev for d in data]

def calculate_dcritical(n, x):
  from math import sqrt, pi, exp
  result = 0.0
  for k in range(1, n):
    a = (2.0 * k - 1)**2
    b = exp(-a * pi**2 / (8*x**2))
    result += b

  result *= sqrt(2.0 * pi) / x
  return result

if __name__ == '__main__':

  #n = 100000
  #alpha = 0.9
  #print 1.0 - calculate_dcritical(n, alpha)

  from random import gauss, uniform
  from matplotlib import pylab

  mean = 20
  sigma = 5
  n = 1000
  data = [uniform(mean, sigma) for i in range(n)]

  data = standardize(data)

  t, f = empirical_distribution_function(data)

  d = distance(t, f)


  from math import sqrt
  D = kolmogorov_smirnov_test(t, f)
  n = len(data)
  Dc = 1.36 / sqrt(n)
  print D, Dc, D < Dc

  pylab.plot(t, f)
  pylab.plot(t, [standard_normal_cdf(tt) for tt in t])
  pylab.show()

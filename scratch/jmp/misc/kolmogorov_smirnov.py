



def standard_normal_cdf(x):
  from math import erf, sqrt
  return (1.0 / 2.0) * (erf(x / sqrt(2.0)) + 1.0)


def distance(t, f):
  f2 = [standard_normal_cdf(tt) for tt in t]
  return [ff2 - ff for ff2, ff in zip(f2, f)]

def supremum(x):
  return max([abs(xx) for xx in x])

def kolmogorov_smirnov_test(t, f):
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

if __name__ == '__main__':

  from random import gauss, uniform
  from matplotlib import pylab

  mean = 10
  sigma = 1
  data = [uniform(mean, sigma) for i in range(1000)]

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

from __future__ import division


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

class Test(object):
  def __init__(self):
    pass

  def run(self):
    for i in range(100, 1000):
      self.tst_normal_for_sample_size(i)
    print 'OK'

    for i in range(500, 1000):
      self.tst_uniform_for_sample_size(i)
    print 'OK'

  def tst_normal_for_sample_size(self, n):
    from dials.algorithms.background import ks_is_normally_distributed
    from random import gauss, uniform
    from dials.array_family import flex

    m = 20
    s = 5
    normal_data = [gauss(m, s) for i in range(n)]

    alpha = 0.05

    assert(ks_is_normally_distributed(flex.double(normal_data), alpha))

  def tst_uniform_for_sample_size(self, n):
    from dials.algorithms.background import ks_is_normally_distributed
    from random import gauss, uniform
    from dials.array_family import flex

    a = 15
    b = 25
    uniform_data = [uniform(a, b) for i in range(n)]
    alpha = 0.2

    #print n
    try:
      assert(not ks_is_normally_distributed(flex.double(uniform_data), alpha))
    except Exception:
      from math import sqrt
      n = len(uniform_data)
      print kolmogorov_smirnov_test(uniform_data), 1.36 / sqrt(n)
      #from matplotlib import pylab

      #uniform_data = standardize(uniform_data)

      #t, f = empirical_distribution_function(uniform_data)

      #pylab.plot(t, f)
      #pylab.plot(t, [standard_normal_cdf(tt) for tt in t])
      #pylab.show()

if __name__ == '__main__':
  test = Test()
  test.run()

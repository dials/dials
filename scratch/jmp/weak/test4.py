
from __future__ import division
from math import log, exp, factorial

def dt1(c, b, s, t1, t2, mu0, mu1, mu2):
  sumb = mu0 * sum(b)
  sums = mu2 * sum(s)
  sum2 = 0
  if (t1 < -20):
    return 0
  for i in range(len(c)):
    den = exp(t1)*(b[i]*mu0-s[i]*mu2) + exp(t2)*(s[i]*mu0-b[i]*mu1)
    try:
      assert(den > 0)
    except Exception, e:
      print "Den:"
      print den
      print t1, t2
      print exp(t1), b[i]*mu0, s[i]*mu2
      print exp(t2), s[i]*mu0, b[i]*mu1
      raise e
    sum2 += c[i]*(b[i]*mu0-s[i]*mu2) / den
  return exp(t1) * (sum2 - sumb + sums)

def dt2(c, b, s, t1, t2, mu0, mu1, mu2):
  sumb = mu1 * sum(b)
  sums = mu0 * sum(s)
  sum2 = 0
  if (t2 < -20):
    return 0
  for i in range(len(c)):
    den = exp(t1)*(b[i]*mu0-s[i]*mu2) + exp(t2)*(s[i]*mu0-b[i]*mu1)
    try:
      assert(den > 0)
    except Exception, e:
      print "Den:"
      print den
      print exp(t1), b[i]*mu0, s[i]*mu2
      print exp(t2), s[i]*mu0, b[i]*mu1
      raise e
    sum2 += c[i]*(s[i]*mu0-b[i]*mu1) / den
  return exp(t2) * (sum2 - sums + sumb)


def iterate(c, b, s, t1, t2, mu0, mu1, mu2):
  d1 = dt1(c, b, s, t1, t2, mu0, mu1, mu2)
  d2 = dt2(c, b, s, t1, t2, mu0, mu1, mu2)
  print d1, d2
  return t1 + d1, t2 + d2

def solve(c, b, s, t1, t2, mu0, mu1, mu2):
  for i in range(10):
    t1_new, t2_new = iterate(c, b, s, t1, t2, mu0, mu1, mu2)
    t1 = t1_new
    t2 = t2_new

    x = mu0*exp(t1) - mu1*exp(t2)
    y = mu0*exp(t2) - mu2*exp(t1)
    print t1, t2, x, y


c = [0, 1, 0, 0, 0, 0, 3, 1, 3, 3, 6, 6, 4, 1, 4, 0, 2, 0, 1, 1]
# c = [0, 1, 0, 0, 0, 0, 0, 0, 2, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
b = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
s = [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0]
b = [bb / sum(b) for bb in b]
s = [ss / sum(s) for ss in s]

k1 = 0
k2 = 9.0 / 20.0
mu0 = 1.0 / (1.0 - k1*k2)
mu1 = k1 / (1.0 - k1*k2)
mu2 = k2 / (1.0 - k1*k2)

solve(c, b, s, 0, 0, mu0, mu1, mu2)

# exit(0)

LL = []
t1 = 0
for i in range(100):
  t2 = -5.0 + i / 10.0
  sumc = sum([log(factorial(cc)) for cc in c])
  sum1 = 0
  sum2 = 0
  for j in range(len(c)):
    bbpss = exp(t1)*(b[j]*mu0-s[j]*mu2) + exp(t2)*(s[j]*mu0-b[j]*mu1)
    sum1 += c[j]*log(bbpss)
    sum2 += bbpss
  LL.append(sum1 - sumc - sum2)
from matplotlib import pylab
pylab.plot(LL)
pylab.show()

exit(0)

from dials.array_family import flex
im1 = flex.double(flex.grid(100, 100))
im2 = flex.double(flex.grid(100, 100))

for j in range(100):
  for i in range(100):
    t1 = -90.0 + i
    t2 = -90.0 + j
    d1 = dt1(c, b, s, t1, t2, mu0, mu1, mu2)
    d2 = dt2(c, b, s, t1, t2, mu0, mu1, mu2)
    im1[j,i] = d1
    im2[j,i] = d2

print flex.min(im1), flex.max(im1)
print flex.min(im2), flex.max(im2)

from matplotlib import pylab
pylab.imshow(im1.as_numpy_array())
pylab.contour(im1.as_numpy_array(), levels=[0])
pylab.contour(im2.as_numpy_array(), levels=[0])
pylab.show()

# print mu0, mu1, mu2

# b = [bb / sum(b) for bb in b]
# s = [ss / sum(s) for ss in s]

# K1 = 0
# K2 = 9.0 / 20.0
# B0 = 1
# S0 = 1
# t2 = 1
# L = []
# for i in range(100):
#   t1 = -98.0 + i
#   sum_c = sum([log(factorial(cc)) for cc in c])
#   sum_1 = 0
#   sum_2 = 0
#   for i in range(len(c)):
#     bbpss = (exp(t1) - K1*S0)*b[i]+(exp(t2) - K2*B0)*s[i]
#     sum_1 += c[i]*log(bbpss)
#     sum_2 += bbpss
#   LL = sum_1 - sum_c - sum_2
#   L.append(LL)


# from matplotlib import pylab
# pylab.plot(L)
# pylab.show()

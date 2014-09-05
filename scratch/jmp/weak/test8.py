
from __future__ import division

c = [0, 1, 0, 0, 0, 0, 3, 1, 3, 3, 6, 6, 4, 1, 4, 0, 2, 0, 1, 1]
c = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
# c = [0, 0, 0, 0, 0, 0, 0, 0, 2, 10, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0]
b = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
s = [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0]
b = [bb / sum(b) for bb in b]
s = [ss / sum(s) for ss in s]

bj = b[9]
sj = s[9]
mu = sum(s) / sum(b)
lb = sum(c) / sum(b)
print 1 - (sj-bj)/((sj - bj*mu) + bj*lb), 9.0 / 20.0, 9 / (20 - 9)
exit(0)


X = []
D = []
sumc = sum(c)
mu = sum(s) / sum(b)
lb = sum(c) / sum(b)
print mu, lb
for j in range(100):
  S = j*(sumc -0) / 100.0
  X.append(S)
  sum1 = 0
  sum2 = 0
  for i in range(len(c)):
    den = (S*(s[i]-b[i]*mu) + b[i]*lb)
    sum1 += c[i]*(s[i] - mu*b[i]) / den
    # den = (S*(s[i] - b[i]) + b[i]*sumc)
    # print den
    # DD += c[i]*(S*(s[i] + b[i]) - s[i]*sumc) / den
  D.append(sum1)

ind = sorted(range(len(D)), key=lambda x: abs(D[x]))[0]
print "Min: ", D[ind], X[ind]

from matplotlib import pylab
pylab.plot(X, D)
pylab.show()



# S = 1

# for j in range(10):

#   f = 0
#   fp = 0
#   fpp = 0
#   for i in range(len(c)):
#     u = s[i]-b[i]*mu
#     v = S*(s[i]-b[i]*mu) + b[i]*lb
#     f += c[i]*u/v
#     fp += -c[i]*u*u/(v*v)
#     fpp += 2*c[i]*u*u*u/(v*v*v)
#   # S = S - 2*f*fp / (2*fp*fp - f*fpp)
#   S = S - f / (fp - (f*fpp/(2*fp)))
#   print S, f, fp, fpp
#   # print S

# from math import log, exp, atan, pi, factorial, tan, sin, cos

# sumb = sum(b)
# sums = sum(s)
# sumc = sum(c)
# sumcf = sum([log(factorial(int(cc))) for cc in c])

# from dials.array_family import flex
# LL = flex.double(flex.grid(100, 100))
# t0 = 0
# t1 = pi/2
# LL2 = []
# LL3 = []
# for j in range(100):
#   # v = -50.0 + j
#   # theta = (t1 - t0)*(atan(v) / pi) + (t1 + t0) / 2.0

#   dt = 0.05 * pi / 2
#   theta = dt + j* (pi / 2)*(1 - 0.1) / 100.0
#   sint = sin(theta)
#   cost = cos(theta)
#   for i in range(100):
#     r = 0.05 + i
#     # u = -10.0 + i / 5.0
#     # r = exp(u)
#     print r
#     sumt = 0
#     for k in range(len(c)):
#       sumt += c[k]*log(b[k]*cost + s[k]*sint)
#     # L = u*sumc + sumt - sumcf - exp(u)*(cost*sumb+sint*sums)
#     L = log(r)*sumc + sumt - sumcf - r*(cost*sumb+sint*sums)
#     LL[j,i] = L
#     if i == 64:
#       LL2.append(L)
#     if j == 50:
#       LL3.append(L)

# from matplotlib import pylab, cm
# pylab.plot(LL2)
# pylab.show()
# pylab.plot(LL3)
# pylab.show()
# pylab.imshow(LL.as_numpy_array(), origin='bottom', cmap=cm.hot)
# pylab.show()

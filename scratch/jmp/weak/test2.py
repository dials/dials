

from math import log, factorial, exp


def f(c, b, s, B, S):
  LL = 0
  for j in range(len(c)):
    LL += c[j] * log(B*b[j]+S*s[j]) - log(factorial(c[j])) - (B*b[j]+S*s[j])
  return LL

def dfdb(c, b, s, B, S):
  LL = 0
  for j in range(len(c)):
    LL += c[j] * b[j] / (B*b[j] + S*s[j]) - b[j]
  return LL

def dfds(c, b, s, B, S):
  LL = 0
  for j in range(len(c)):
    LL += c[j] * s[j] / (B*b[j] + S*s[j]) - s[j]
  return LL

def d2fdb2(c, b, s, B, S):
  LL = 0
  for j in range(len(c)):
    LL += -c[j] * b[j] * b[j] / (B*b[j] + S*s[j])**2
  return LL

def d2fds2(c, b, s, B, S):
  LL = 0
  for j in range(len(c)):
    LL += -c[j] * s[j] * s[j] / (B*b[j] + S*s[j])**2
  return LL

def d2fdbds(c, b, s, B, S):
  LL = 0
  for j in range(len(c)):
    LL += -c[j] * b[j] * s[j] / (B*b[j] + S*s[j])**2
  return LL

c = [0, 1, 0, 0, 0, 0, 3, 1, 3, 3, 6, 6, 4, 1, 4, 0, 2, 0, 1, 1]
b = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
s = [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0]

S = 0

LL = []
for i in range(100):
  B = 0.1 + i / 10.0
  LLL = 0
  for j in range(len(c)):
    LLL += c[j] * log(B*b[j]+S*s[j]) - log(factorial(c[j])) - (B*b[j]+S*s[j])
  LL.append(LLL)

BB = []
LLA = []
B0 = 1
A1 = sum([log(factorial(cc)) for cc in c])
A2 = sum(b)
A3 = sum(s)
K1 = -d2fdb2(c, b, s, B0, S) * B0**2 / A2
K2 = (1.0 / B0)*exp(A1 + A2*B0 + f(c, b, s, B0, S) / (d2fdb2(c, b, s, B0, S)*B0**2))
K3 = 1.0
for i in range(100):
  B = 0.1 + i / 10.0
  LLL = K1 * log(B*K2 + S*K3) - A1 - A2*B - A3*S
  LLA.append(LLL)
  BB.append(B)

Q = []
F = f(c, b, s, B0, S)
DF = dfdb(c, b, s, B0, S)
D2F = d2fdb2(c, b, s, B0, S)
K1 = D2F
K2 = B0 - DF / D2F
K3 = F - DF*DF / (2*D2F)
for i in range(100):
  B = 0.1 + i / 10.0
  Q.append(K3 + (K1 / 2) *(B - K2)**2)

from matplotlib import pylab
pylab.plot(BB, LL)
pylab.plot(BB, LLA)
# pylab.plot(BB, Q)
pylab.show()


def db(c, b, s, B, S):
  sumb = sum(b)
  sumc = 0
  for i in range(len(c)):
    sumc += c[i]*b[i] / (B*b[i] + S*s[i])
  return sumc - sumb

def ds(c, b, s, B, S):
  sums = sum(s)
  sumc = 0
  for i in range(len(c)):
    sumc += c[i]*s[i] / (B*b[i] + S*s[i])
  return sumc - sums

def d2db2(c, b, s, B, S):
  sumc = 0
  for i in range(len(c)):
    sumc += c[i]*b[i]*b[i] / (B*b[i]+S*s[i])**2
  return -sumc

def d2ds2(c, b, s, B, S):
  sumc = 0
  for i in range(len(c)):
    sumc += c[i]*s[i]*s[i] / (B*b[i]+S*s[i])**2
  return -sumc

def d2dbds(c, b, s, B, S):
  sumc = 0
  for i in range(len(c)):
    sumc += c[i]*s[i]*b[i] / (B*b[i]+S*s[i])**2
  return -sumc

def d3db3(c, b, s, B, S):
  sumc = 0
  for i in range(len(c)):
    sumc += c[i]*b[i]*b[i]*b[i] / (B*b[i]+S*s[i])**3
  return sumc

def d3db2ds(c, b, s, B, S):
  sumc = 0
  for i in range(len(c)):
    sumc += c[i]*b[i]*b[i]*s[i] / (B*b[i]+S*s[i])**3
  return sumc

def d3dbds2(c, b, s, B, S):
  sumc = 0
  for i in range(len(c)):
    sumc += c[i]*b[i]*s[i]*s[i] / (B*b[i]+S*s[i])**3
  return sumc

def d3ds3(c, b, s, B, S):
  sumc = 0
  for i in range(len(c)):
    sumc += c[i]*s[i]*s[i]*s[i] / (B*b[i]+S*s[i])**3
  return sumc

c = [0, 1, 0, 0, 0, 0, 3, 1, 3, 3, 6, 6, 4, 1, 4, 0, 2, 0, 1, 1]
c = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
b = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
s = [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0]
b = [bb / sum(b) for bb in b]
s = [ss / sum(s) for ss in s]

from scipy.optimize import minimize

def f(x, c, b, s):
  from math import log, factorial
  print x
  sum_c = sum([log(factorial(cc)) for cc in c])
  sum_b = sum(b)
  sum_s = sum(s)
  sum_l = 0
  for i in range(len(c)):
    assert(x[0]*b[i]+x[1]*s[i] > 0)
    sum_l += c[i]*log(x[0]*b[i]+x[1]*s[i])
  return -(sum_l - sum_c - x[0]*sum_b - x[1]*sum_s)

def J(x, c, b, s):
  import numpy
  return numpy.array([-db(c, b, s, x[0], x[1]), -ds(c, b, s, x[0], x[1])])

def C1(x):
  return x[0]

def C2(x):
  return x[0] + 0.45*x[1]

def C1J(x):
  import numpy
  return numpy.array([1, 0])

def C2J(x):
  import numpy
  return numpy.array([1, 0.45])

con = [{
  'type' : 'ineq',
  'fun' : C1,
  'jac' : C1J
}, {
  'type' : 'ineq',
  'fun' : C2,
  'jac' : C2J
}]

x0 = [1, 1]
print minimize(f, x0, args=(c, b, s), method="SLSQP", constraints=con, jac=J)
# from scitbx import matrix

# B = 1
# S = 1

# for i in range(10):

#   F = matrix.col((
#     db(c,b,s,B,S), 
#     ds(c,b,s,B,S)))

#   FP = matrix.sqr((
#     d2db2(c,b,s,B,S), d2dbds(c,b,s,B,S),
#     d2dbds(c,b,s,B,S), d2ds2(c,b,s,B,S)))


#   AV = -FP.inverse() * F
#   BV = FP.inverse() * FPP * AV*AV

#   # X = matrix.col((B,S))
#   # X = X - H.inverse() * J
#   # B = X[0]
#   # S = X[1]
#   print X

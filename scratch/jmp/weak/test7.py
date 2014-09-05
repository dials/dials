from __future__ import division
import numpy as np
from scipy import optimize
import pylab as pl
from mpl_toolkits.mplot3d.axes3d import Axes3D

# c = [0, 1, 0, 0, 0, 0, 3, 1, 3, 3, 6, 6, 4, 1, 4, 0, 2, 0, 1, 1]
c = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
b = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
s = [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0]
b = [bb / sum(b) for bb in b]
s = [ss / sum(s) for ss in s]

# H = np.array([  [2.,0.],
#         [0.,8.] ])
# c = np.array(   [0,-32])
# c0 = 64

A = np.array([
  [1.0, 0.0],
  [0.45, 1.0]])

# A = np.array([  [ 1.,1.],
#         [-1.,2.],
#         [-1.,0.],
#         [0.,-1.],
#         [0., 1.] ])

C = np.array([0, 0])
# b = np.array(   [7.,4.,0.,0.,4.] )

x0 = np.array([1e-7, 1e-7])

def func(B, S):
  from math import log, factorial
  sumc = sum([log(factorial(cc)) for cc in c])
  sumb = sum(b)
  sums = sum(s)
  suml = 0
  for i in range(len(c)):
    assert(B*b[i] + S*s[i] > 1e-10)
    suml += c[i] * log(B*b[i] + S*s[i])
  return -(suml - sumc - B*sumb - S*sums)

def objective(x0):
  print 1, x0
  return func(x0[0], x0[1])

def db(B, S):
  sumb = sum(b)
  sumc = 0
  for i in range(len(c)):
    assert(B*b[i] + S*s[i] > 1e-10)
    sumc += c[i]*b[i] / (B*b[i] + S*s[i])
  return sumc - sumb

def ds(B, S):
  sums = sum(s)
  sumc = 0
  for i in range(len(c)):
    assert(B*b[i] + S*s[i] > 1e-10)
    sumc += c[i]*s[i] / (B*b[i] + S*s[i])
  return sumc - sums

def jacobian(x0):
  print 2, x0
  import numpy
  return numpy.array((
    -db(x0[0], x0[1]),
    -ds(x0[0], x0[1])))

# def objective(x,sign=1.):
#     return sign*(0.5*np.dot(x.T,np.dot(H,x))+ np.dot(c,x) + c0)

# def jacobian(x,sign=1.):
#     return sign*(np.dot(x.T,H) + c)

def constraints(x):
  CC = np.dot(A,x) - C
  print "C: ", x, CC
  return CC

cons = (
  {
    'type':'ineq',
    'fun':constraints
  })

def solve():
    res_cons = optimize.minimize(
      objective,
      x0,
      jac=jacobian,
      constraints=cons,
      method='SLSQP',
      options={'disp':False})

    print '\nConstrained:'
    print res_cons

    # print '\nUnconstrained:'
    # print res_uncons

    # x1,x2 = res_cons['x']
    # f = res_cons['fun']

    # x1_unc,x2_unc = res_uncons['x']
    # f_unc = res_uncons['fun']

    # xgrid = np.mgrid[-2:4:0.1,1.5:5.5:0.1]
    # xvec = xgrid.reshape(2,-1).T
    # F = np.vstack([objective(xi) for xi in xvec]).reshape(xgrid.shape[1:])

    # ax = pl.axes(projection='3d')
    # ax.hold(True)
    # ax.plot_surface(xgrid[0],xgrid[1],F,rstride=1,cstride=1,cmap=pl.cm.jet,
    #     shade=True,alpha=0.9,linewidth=0)
    # ax.plot3D(np.atleast_1d(x1),np.atleast_1d(x2),np.atleast_1d(f),
    #     'og',mec='w',label='Constrained minimum')
    # ax.plot3D(np.atleast_1d(x1_unc),np.atleast_1d(x2_unc),np.atleast_1d(f_unc),
    #     'oy',mec='w',label='Unconstrained minimum')

    # ax.legend(fancybox=True,numpoints=1)
    # ax.set_xlabel('x1')
    # ax.set_ylabel('x2')
    # ax.set_zlabel('F')
    # pl.show()


solve()


import numpy as np
from scipy.stats import nbinom

def _ll_nb2(y, X, beta, alph):
  mu = np.exp(np.dot(X, beta))
  size = 1 / alph
  prob = size / (size + mu)
  ll = nbinom.logpmf(y, size, prob)
  return ll

from statsmodels.base.model import GenericLikelihoodModel

class NBin(GenericLikelihoodModel):
     def __init__(self, endog, exog, **kwds):
         super(NBin, self).__init__(endog, exog, **kwds)
     def nloglikeobs(self, params):
         alph = params[-1]
         beta = params[:-1]
         ll = _ll_nb2(self.endog, self.exog, beta, alph)
         return -ll
     def fit(self, start_params=None, maxiter=10000, maxfun=5000, **kwds):
         if start_params == None:
             # Reasonable starting values
             start_params = np.append(np.zeros(self.exog.shape[1]), .5)
             start_params[0] = np.log(self.endog.mean())
         return super(NBin, self).fit(start_params=start_params,
                                      maxiter=maxiter, maxfun=maxfun,
                                      **kwds)


def nb_fit(y):
  import numpy

  y = numpy.array([[yy] for yy in y])
  X = numpy.array([[1.0] for yy in y])
  mod = NBin(y, X)
  res = mod.fit()
  return tuple(res.params)

from numpy.random import poisson, negative_binomial

y = list(poisson(100, 100))
y = list(negative_binomial(1, 0.9, 100))
from matplotlib import pylab
pylab.hist(y)
pylab.show()
print y
print nb_fit(y)

# y = numpy.array(
# [[ 4.0],
#  [ 9.0],
#  [ 3.0],
#  [ 9.0],
#  [ 1.0]])

# X = numpy.array(
# [[ 1.0],
#  [ 1.0],
#  [ 1.0],
#  [ 1.0],
#  [ 1.0]])
# print y[:5]

# print X[:5]


# mod = NBin(y, X)
# res = mod.fit()

# print res.params

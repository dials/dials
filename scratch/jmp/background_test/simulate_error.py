from __future__ import division

def generate_distribution(mean, sigma, npixels):
  from random import randint
  from scitbx.array_family import flex
  values = [0 for i in range(npixels)]
  for j in range(int(mean * npixels)):
    x = randint(0, npixels-1)
    values[x] += 1
  return flex.double(values)

def index_of_dispersion(data):
  from scitbx.array_family import flex
  mv = flex.mean_and_variance(data)
  mean = mv.mean()
  var = mv.unweighted_sample_variance()
  return var / mean

def index_of_dispersion_expected_variance(npixels):
  return 2.0 / (npixels - 1)

def is_poisson_distributed(data):
  from math import sqrt
  return index_of_dispersion(data) < 1.0 + \
    sqrt(index_of_dispersion_expected_variance(len(data)))

if __name__ == '__main__':
  from math import sqrt
  from dials.algorithms.background import normal_expected_n_sigma, is_normally_distributed, maximum_n_sigma
  import numpy
  from random import sample
  from scitbx.array_family import flex

  mean = 10.0
  sigma = 5
  npixels = 1000
  values = generate_distribution(mean, sigma, npixels)

  ind = sample(range(1000), 800)
  values = [float(values[i]) for i in ind]

  cmean = sum(values) / len(values)
  cvar = sum([(v - cmean)**2 for v in values]) / (len(values) - 1)
  csdev = sqrt(cvar)
  cerr = csdev / sqrt(len(values))

  print "Mean: %f, Error: %f" % (cmean, cerr)


#  iod = []
#  for k in range(10000):
#    values = numpy.random.poisson(10, 1000)
#    mean = numpy.mean(values)
#    var = numpy.var(values)
#    i = var / mean
#    iod.append(i)
#
#

##  iod = []
##  for i in range(10000):

##    mean = 10.0
##    sigma = 5
##    npixels = 1000
##    values = generate_distribution(mean, sigma, npixels)

##    iod.append(index_of_dispersion(values))

#

#  from matplotlib import pylab
##  low = 1.0 - sqrt(index_of_dispersion_expected_variance(1000))
##  high = 1.0 + sqrt(index_of_dispersion_expected_variance(1000))
##  pylab.plot(sorted(iod))
##  pylab.axhline(y=low)
##  pylab.axhline(y=high)
##  pylab.show()
#  pylab.hist(iod)
#  pylab.show()
##  high = sqrt(index_of_dispersion_expected_variance(1000))
##  pylab.plot(sorted([abs(1.0 - i) for i in iod]))
##  pylab.axhline(y=high)
##  pylab.show()


#  #print "Is Normal: ", maximum_n_sigma(flex.double(iod)), normal_expected_n_sigma(len(iod)), is_normally_distributed(flex.double(iod))
##  print "Is Poisson: ", is_poisson_distributed(values)

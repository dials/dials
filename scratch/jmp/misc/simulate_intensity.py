from __future__ import division

def simulate(n, size):
  from scitbx.array_family import flex
  from scitbx.random import variate, poisson_distribution
  shoeboxes = []

  B = 10

  # Generate shoeboxes with uniform random background
  for l in range(n):
    sbox = flex.double(flex.grid(size),0)
    g = variate(poisson_distribution(mean = B))
    for k in range(size[0]):
      for j in range(size[1]):
        for i in range(size[2]):
          sbox[k, j, i] += g.next()
    shoeboxes.append(sbox)

  # Calculate the Intensity (should be zero)
  import random
  I_cal = []
  mean = []
  for i in range(len(shoeboxes)):
    nn = len(shoeboxes[i])
    mm = int(1.0 * nn)
    index = flex.size_t(random.sample(range(nn), mm))
    assert(len(set(index)) == mm)
    data = shoeboxes[i].select(index)
    II = flex.sum(data)
    #II = flex.mean(data)
    BB = mm * B
    #BB = B
    I_cal.append(II - BB)
    m = flex.mean(data)
    mean.append(m)
  I_cal = flex.double(I_cal)

  mv = flex.mean_and_variance(flex.double(mean))
  print mv.mean() - B, mv.unweighted_sample_variance()
  v1 = B / (size[0] * size[1] * size[2])
  v2= B * (size[0] * size[1] * size[2])
  print v1
  print v2
  print I_cal[0]

  from math import sqrt
  Z = (I_cal - 0) / sqrt(v2)


  # Return the mean and standard deviation
  mv = flex.mean_and_variance(Z)
  return mv.mean(), mv.unweighted_sample_variance()

if __name__ == '__main__':
  from scitbx.array_family import flex
  mean, sdev = [], []
  for i in range(1):
    m, s = simulate(1000, (10, 10, 10))
    print "Mean: %f, Var: %f" % (m, s)
    mean.append(m)
    sdev.append(s)
#  mv = flex.mean_and_variance(flex.double(mean))
  #print "Mean of Means: ", mv.mean()
  #print "SDev of Means: ", mv.unweighted_sample_standard_deviation()

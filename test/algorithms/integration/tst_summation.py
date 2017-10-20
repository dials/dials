

from __future__ import absolute_import, division

class Test(object):

  def __init__(self):
    from os.path import join, isfile
    import libtbx.load_env
    try:
      dials_regression = libtbx.env.dist_path('dials_regression')
    except KeyError:
      print 'SKIP: dials_regression not configured'
      exit(0)

    # The base path
    path = join(dials_regression, 'integration_test_data', 'simulated')

    # The paths to the reflection files
    self.refl_filenames = [
      join(path, 'simulated_n10000_i0_b10.pickle'),
      join(path, 'simulated_n10000_i0_b100.pickle'),
      join(path, 'simulated_n10000_i0_b1000.pickle'),
      join(path, 'simulated_r_n10000_i0_b1000.pickle'),

      join(path, 'simulated_n10000_i10_b0.pickle'),
      join(path, 'simulated_n10000_i100_b0.pickle'),
      join(path, 'simulated_n10000_i1000_b0.pickle'),
      join(path, 'simulated_n10000_i10000_b0.pickle'),
      join(path, 'simulated_r_n10000_i10000_b0.pickle'),

      join(path, 'simulated_n10000_i10_b10.pickle'),
      join(path, 'simulated_n10000_i100_b10.pickle'),
      join(path, 'simulated_n10000_i1000_b10.pickle'),
      join(path, 'simulated_n10000_i10000_b10.pickle'),
      join(path, 'simulated_r_n10000_i10000_b10.pickle'),

      join(path, 'simulated_n10000_i10_b100.pickle'),
      join(path, 'simulated_n10000_i100_b100.pickle'),
      join(path, 'simulated_n10000_i1000_b100.pickle'),
      join(path, 'simulated_n10000_i10000_b100.pickle'),
      join(path, 'simulated_r_n10000_i10000_b100.pickle'),
    ]

    # Check the files exist
    for filename in self.refl_filenames:
      if not isfile(filename):
        print 'SKIP: simulated test data does not exist'
        print 'Generate by running the following commands:'
        print ' cd dials_regression/integration_test_data/simulated'
        print ' ./simulate'
        exit(0)

  def run(self):
    from dials.array_family import flex
    for filename in self.refl_filenames:
      refl = flex.reflection_table.from_pickle(filename)
      self.test_for_reflections(refl)

  def mv3n_tolerance_interval(self, x):
      from math import sqrt, pi, erf, exp
      g32 = sqrt(pi) / 2.0
      x2 = x / 2.0
      srx2 = sqrt(x2)
      return (g32 * erf(srx2) - srx2 * exp(-x2)) / g32

  def test_for_reflections(self, refl):
    from dials.algorithms.integration.sum import IntegrationAlgorithm
    from dials.array_family import flex
    from dials.algorithms.statistics import \
      kolmogorov_smirnov_test_standard_normal

    # Get the calculated background and simulated background
    B_sim = refl['background.sim.a'].as_double()
    I_sim = refl['intensity.sim'].as_double()
    I_exp = refl['intensity.exp']

    # Set the background as simulated
    shoebox = refl['shoebox']
    for i in range(len(shoebox)):
      bg = shoebox[i].background
      ms = shoebox[i].mask
      for j in range(len(bg)):
        bg[j] = B_sim[i]


    # Integrate
    integration = IntegrationAlgorithm()
    integration(refl)
    I_cal = refl['intensity.sum.value']
    I_var = refl['intensity.sum.variance']

    # Only select variances greater than zero
    mask = I_var > 0
    I_cal = I_cal.select(mask)
    I_var = I_var.select(mask)
    I_sim = I_sim.select(mask)
    I_exp = I_exp.select(mask)

    # Calculate the z score
    perc = self.mv3n_tolerance_interval(3*3)
    Z = (I_cal - I_exp) / flex.sqrt(I_var)
    mv = flex.mean_and_variance(Z)
    Z_mean = mv.mean()
    Z_var = mv.unweighted_sample_variance()
    print "Z: mean: %f, var: %f" % (Z_mean, Z_var)

    # Do the kolmogorov smirnov test
    D, p  = kolmogorov_smirnov_test_standard_normal(Z)
    print "KS: D: %f, p-value: %f" % (D, p)

    # FIXME Z score should be a standard normal distribution. When background is
    # the main component, we do indeed see that the z score is in a standard
    # normal distribution. When the intensity dominates, the variance of the Z
    # scores decreases indicating that for increasing intensity of the signal,
    # the variance is over estimated.
    assert(abs(Z_mean) <= 3 * Z_var)


    #from matplotlib import pylab
    #pylab.hist(Z, 20)
    #pylab.show()

    #Z_I = sorted(Z)
    ##n = int(0.05 * len(Z_I))
    ##Z_I = Z_I[n:-n]
    ##mv = flex.mean_and_variance(flex.double(Z_I))
    ##print "Mean: %f, Sdev: %f" % (mv.mean(), mv.unweighted_sample_standard_deviation())
    #edf = [float(i+1) / len(Z_I) for i in range(len(Z_I))]
    #cdf = [0.5 * (1.0 + erf(z / sqrt(2.0))) for z in Z_I]

    print 'OK'

if __name__ == '__main__':

  test = Test()
  test.run()

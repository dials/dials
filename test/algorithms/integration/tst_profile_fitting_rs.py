

from __future__ import division

class Test(object):

  def __init__(self):
    from os.path import join, isfile
    import libtbx.load_env
    try:
      dials_regression = libtbx.env.dist_path('dials_regression')
    except KeyError, e:
      print 'SKIP: dials_regression not configured'
      exit(0)

    # The base path
    path = join(dials_regression, 'integration_test_data', 'simulated')

    # Set the experiment filename
    expr_filename = join(path, 'experiments.json')

    # The reference spots filesname
    reference_filename = join(path, 'simulated_n10000_i10000_b0.pickle')

    # The paths to the reflection files
    self.refl_filenames = [
      #join(path, 'simulated_n10000_i0_b10.pickle'),
      #join(path, 'simulated_n10000_i0_b100.pickle'),
      #join(path, 'simulated_n10000_i0_b1000.pickle'),
      #join(path, 'simulated_r_n10000_i0_b1000.pickle'),

      #join(path, 'simulated_n10000_i10_b0.pickle'),
      #join(path, 'simulated_n10000_i100_b0.pickle'),
      #join(path, 'simulated_n10000_i1000_b0.pickle'),
      #join(path, 'simulated_r_n10000_i10000_b0.pickle'),

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

    # Load the experiments
    from dxtbx.model.experiment.experiment_list import ExperimentListFactory
    experiments = ExperimentListFactory.from_json_file(expr_filename,
                                                       check_format=False)
    assert(len(experiments) == 1)
    self.experiments = experiments

    from dials.algorithms.profile_model.model_list import ProfileModelList
    from dials.algorithms.profile_model.gaussian_rs import ProfileModel
    self.profile_model = ProfileModelList()
    self.profile_model.append(
      ProfileModel(
        sigma_b=0.024,
        sigma_m=0.044,
        n_sigma=3,
        deg=True))

    # Load the reference spots
    from dials.array_family import flex
    self.reference = flex.reflection_table.from_pickle(reference_filename)
    mask = flex.bool(len(self.reference), True)
    value = self.reference.flags.reference_spot
    self.reference.set_flags(mask, value)

  def run(self):
    from dials.array_family import flex
    #self.test_for_reference()
    for filename in self.refl_filenames:
      refl = flex.reflection_table.from_pickle(filename)
      refl['id'] = flex.int(len(refl),0)
      self.test_for_reflections(refl, filename)

  def mv3n_tolerance_interval(self, x):
      from math import sqrt, pi, erf, exp
      g32 = sqrt(pi) / 2.0
      x2 = x / 2.0
      srx2 = sqrt(x2)
      return (g32 * erf(srx2) - srx2 * exp(-x2)) / g32

  def test_for_reference(self):
    from dials.algorithms.integration.fitrs import IntegrationAlgorithm
    from dials.array_family import flex
    from math import sqrt, pi

    # Integrate
    integration = IntegrationAlgorithm(
      grid_size=4,
      threshold=0.02,
      frame_interval=100,
      n_sigma=5,
      sigma_b=0.024 * pi / 180.0,
      sigma_m=0.044 * pi / 180.0
    )

    # Integrate the reference profiles
    integration(self.experiments, self.reference)
    locator = integration.learner.locate()
    # Check the reference profiles and spots are ok
    #self.check_profiles(integration.learner)

    # Make sure background is zero
    profiles = self.reference['rs_shoebox']
    eps = 1e-7
    for p in profiles:
      assert(abs(flex.sum(p.background) - 0) < eps)
    print 'OK'

    # Only select variances greater than zero
    mask = self.reference.get_flags(self.reference.flags.integrated, all=False)
    assert(mask.count(True) > 0)
    I_cal = self.reference['intensity.prf.value']
    I_var = self.reference['intensity.prf.variance']
    B_sim = self.reference['background.sim.a'].as_double()
    I_sim = self.reference['intensity.sim'].as_double()
    I_exp = self.reference['intensity.exp']
    P_cor = self.reference['profile.correlation']
    X_pos, Y_pos, Z_pos = self.reference['xyzcal.px'].parts()
    I_cal = I_cal.select(mask)
    I_var = I_var.select(mask)
    I_sim = I_sim.select(mask)
    I_exp = I_exp.select(mask)
    P_cor = P_cor.select(mask)

    max_ind = flex.max_index(flex.abs(I_cal-I_sim))
    max_I = I_cal[max_ind]
    max_P = self.reference[max_ind]['rs_shoebox'].data
    max_C = self.reference[max_ind]['xyzcal.px']
    max_S = self.reference[max_ind]['shoebox'].data

    ref_ind = locator.index(max_C)
    ref_P = locator.profile(ref_ind)
    ref_C = locator.coord(ref_ind)

    #def f(I):
      #mask = flex.bool(flex.grid(9,9,9), False)
      #for k in range(9):
        #for j in range(9):
          #for i in range(9):
            #dx = 5 * (i - 4.5) / 4.5
            #dy = 5 * (j - 4.5) / 4.5
            #dz = 5 * (k - 4.5) / 4.5
            #dd = sqrt(dx**2 + dy**2 + dz**2)
            #if dd <= 3:
              #mask[k,j,i] = True

      #mask = mask.as_1d() & (ref_P.as_1d() > 0)
      #p = ref_P.as_1d().select(mask)
      #c = max_P.as_1d().select(mask)
      #return flex.sum((c - I * p)**2 / (I * p))

    #ff = []
    #for I in range(9500, 11500):
      #ff.append(f(I))
    #print 'Old I: ', sorted(range(len(ff)), key=lambda x: ff[x])[0] + 9500

    #from matplotlib import pylab
    #pylab.plot(range(9500, 11500), ff)
    #pylab.show()

    #def estI(I):
      #mask = flex.bool(flex.grid(9,9,9), False)
      #for k in range(9):
        #for j in range(9):
          #for i in range(9):
            #dx = 5 * (i - 4.5) / 4.5
            #dy = 5 * (j - 4.5) / 4.5
            #dz = 5 * (k - 4.5) / 4.5
            #dd = sqrt(dx**2 + dy**2 + dz**2)
            #if dd <= 3:
              #mask[k,j,i] = True

      #mask = mask.as_1d() & (ref_P.as_1d() > 0)
      #p = ref_P.as_1d().select(mask)
      #c = max_P.as_1d().select(mask)
      #v = I * p
      #return flex.sum(c * p / v) / flex.sum(p*p/v)

    #def iterI(I0):
      #I = estI(I0)
      #print I
      #if abs(I - I0) < 1e-3:
        #return I
      #return iterI(I)

    #newI = iterI(10703)#flex.sum(max_P))
    #print "New I: ", newI

    # Calculate the z score
    perc = self.mv3n_tolerance_interval(3*3)
    Z = (I_cal - I_sim) / flex.sqrt(I_var)
    mv = flex.mean_and_variance(Z)
    Z_mean = mv.mean()
    Z_var = mv.unweighted_sample_variance()
    print "Z: mean: %f, var: %f, sig: %f" % (Z_mean, Z_var, sqrt(Z_var))

    # from matplotlib import pylab
    # from mpl_toolkits.mplot3d import Axes3D
    #fig = pylab.figure()
    #ax = fig.add_subplot(111, projection='3d')
    #ax.scatter(X_pos, Y_pos, P_cor)

    #pylab.scatter(X_pos, P_cor)
    #pylab.scatter(Y_pos, P_cor)
    #pylab.scatter(Z_pos, P_cor)
    #pylab.hist(P_cor,100)
    #pylab.scatter(P_cor, (I_cal - I_exp) / I_exp)
    #pylab.hist(Z, 100)
    #pylab.hist(I_cal,100)
    #pylab.hist(I_cal - I_sim, 100)
    #pylab.show()

  def test_for_reflections(self, refl, filename):
    from dials.algorithms.integration.fitrs import IntegrationAlgorithm
    from dials.array_family import flex
    from dials.algorithms.statistics import \
      kolmogorov_smirnov_test_standard_normal
    from os.path import basename
    print basename(filename)

    #refl = self.reference

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
    integration = IntegrationAlgorithm(
      self.experiments,
      self.profile_model,
      grid_size=4,
      threshold=0.00,
    )

    old_size = len(refl)
    refl.extend(self.reference)
    ref_profiles = integration(refl)
    reference = refl[-len(self.reference):]
    refl = refl[:len(self.reference)]
    assert(len(refl) == old_size)
    assert(len(reference) == len(self.reference))
    I_cal = refl['intensity.prf.value']
    I_var = refl['intensity.prf.variance']

    # Check the reference profiles and spots are ok
    #self.check_profiles(integration.learner)
    #self.check_reference(reference)

    #np = reference.locate().size()
    #for i in range(np):
      #profile = reference.locate().profile(i)
      #print "Profile: %d" % i
      #p = (profile.as_numpy_array() * 1000)
      #import numpy as np
      #p = p.astype(np.int)
      #print p

    # Only select variances greater than zero
    mask = refl.get_flags(refl.flags.integrated, all=False)
    assert(mask.count(True) > 0)
    I_cal = I_cal.select(mask)
    I_var = I_var.select(mask)
    I_sim = I_sim.select(mask)
    I_exp = I_exp.select(mask)

    mask = I_var > 0
    I_cal = I_cal.select(mask)
    I_var = I_var.select(mask)
    I_sim = I_sim.select(mask)
    I_exp = I_exp.select(mask)

    # Calculate the z score
    perc = self.mv3n_tolerance_interval(3*3)
    Z = (I_cal - I_sim) / flex.sqrt(I_var)
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
    #assert(abs(Z_mean) <= 3 * Z_var)


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

  def check_profiles(self, learner):
    ''' Check the reference profiles. '''
    from dials.array_family import flex
    from dials.algorithms.image.centroid import centroid_image
    from math import sqrt

    # Get the reference locator
    locator = learner.locate()
    np = locator.size()
    assert(np == 9)
    assert(flex.sum(learner.counts()) == 10000)

    #profile = locator.profile(0)
    #pmax = flex.max(profile)
    #profile = 100 * profile / pmax
    #profile = profile.as_numpy_array()
    #import numpy
    #profile = profile.astype(numpy.int)
    #print profile

    # Check all the profiles
    eps = 1e-7
    n_sigma = 3
    n_sigma2 = 5
    grid_size = 4
    step_size = n_sigma2 / (grid_size + 0.5)
    for i in range(np):
      profile = locator.profile(i)
      assert(abs(flex.sum(profile) - 1.0) < eps)
      centroid = centroid_image(profile)
      m = centroid.mean()
      v = centroid.variance()
      s1 = tuple(sqrt(vv) for vv in v)
      s2 = tuple(ss * step_size for ss in s1)
      assert(all(abs(mm - (grid_size + 0.5)) < 0.25 for mm in m))
      assert(all(abs(ss2 - n_sigma / n_sigma2) < 0.25 for ss2 in s2))

    # Check all the profiles have good correlation coefficients
    cor = locator.correlations()
    assert(all(cor > 0.99))

    # Test passed
    print 'OK'

  def check_reference(self, reference):
    ''' Check the reference spots. '''
    from dials.array_family import flex
    from dials.algorithms.image.centroid import centroid_image
    from math import sqrt

    # Get a load of stuff
    I_sim = reference['intensity.sim']
    I_exp = reference['intensity.exp']
    I_cal = reference['intensity.prf.value']
    I_var = reference['intensity.prf.variance']

    # Get the transformed shoeboxes
    profiles = reference['rs_shoebox']
    n_sigma = 3
    n_sigma2 = 5
    grid_size = 4
    step_size = n_sigma2 / (grid_size + 0.5)
    eps = 1e-7
    for i in range(len(profiles)):
      data = profiles[i].data
      #dmax = flex.max(data)
      #data = 100 * data / dmax
      #p = data.as_numpy_array()
      #p = p.astype(numpy.int)
      #print p
      print flex.sum(data), I_exp[i], I_cal[i]
      #assert(abs(flex.sum(data) - I_exp[i]) < eps)
      centroid = centroid_image(data)
      m = centroid.mean()
      v = centroid.variance()
      s1 = tuple(sqrt(vv) for vv in v)
      s2 = tuple(ss * step_size for ss in s1)
      assert(all(abs(mm - (grid_size + 0.5)) < 0.25 for mm in m))
      assert(all(abs(ss2 - n_sigma / n_sigma2) < 0.25 for ss2 in s2))

    # Calculate Z
    Z = (I_cal - I_exp) / flex.sqrt(I_var)
    mv = flex.mean_and_variance(Z)
    Z_mean = mv.mean()
    Z_var = mv.unweighted_sample_variance()
    print "Z: mean: %f, var: %f" % (Z_mean, Z_var)

    from matplotlib import pylab
    pylab.hist((I_cal - I_exp) / I_exp)
    pylab.show()

if __name__ == '__main__':

  test = Test()
  test.run()

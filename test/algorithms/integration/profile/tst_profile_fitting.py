from __future__ import absolute_import, division
from dials.test.algorithms.integration.profile.tst_profile_helpers import gaussian

def add_poisson_noise(x):
  from dials.array_family import flex
  from numpy.random import poisson
  s = x.accessor()
  y = flex.double(poisson(xx) for xx in x)
  y.reshape(s)
  return y

class Test(object):

  def __init__(self):
    pass

  def run(self):

    self.tst_zero()
    self.tst_identical()
    self.tst_with_no_background()
    self.tst_with_flat_background()
    self.tst_identical_partial()
    self.tst_with_no_background_partial()
    self.tst_with_flat_background_partial()

    #self.tst_deconvolve_zero()

    self.tst_deconvolve_3_with_no_background()
    self.tst_deconvolve_3_with_flat_background()

    self.tst_deconvolve_7_with_no_background()
    self.tst_deconvolve_7_with_flat_background()

  def tst_zero(self):

    from dials.algorithms.integration.fit import ProfileFitter
    from scitbx.array_family import flex
    from numpy.random import seed
    seed(0)

    # Create profile
    p = gaussian((9, 9, 9), 1, (4, 4, 4), (2, 2, 2))
    s = flex.sum(p)
    p = p / s

    # Copy profile
    c = flex.double(flex.grid(9,9,9), 0)
    b = flex.double(flex.grid(9, 9, 9), 0)
    m = flex.bool(flex.grid(9,9,9), True)

    # Fit
    fit = ProfileFitter(c, b, m, p)
    I = fit.intensity()
    V = fit.variance()
    assert fit.niter() < fit.maxiter()

    # Test intensity is the same
    eps = 1e-3
    assert(abs(I[0] - flex.sum(c)) < eps)
    assert(abs(V[0] - I[0]) < eps)

    print 'OK'

  def tst_identical(self):

    from dials.algorithms.integration.fit import ProfileFitter
    from scitbx.array_family import flex
    from numpy.random import seed
    seed(0)

    # Create profile
    p = gaussian((9, 9, 9), 1, (4, 4, 4), (2, 2, 2))
    s = flex.sum(p)
    p = p / s

    # Copy profile
    c = p.deep_copy()
    b = flex.double(flex.grid(9, 9, 9), 0)
    m = flex.bool(flex.grid(9,9,9), True)

    # Fit
    fit = ProfileFitter(c, b, m, p)
    I = fit.intensity()
    V = fit.variance()
    assert fit.niter() < fit.maxiter()

    # Test intensity is the same
    eps = 1e-3
    assert(abs(I[0] - flex.sum(c)) < eps)
    assert(abs(V[0] - I[0]) < eps)

    print 'OK'

  def tst_with_no_background(self):

    from dials.algorithms.integration.fit import ProfileFitter
    from scitbx.array_family import flex
    from numpy.random import seed
    seed(0)

    # Create profile
    p = gaussian((9, 9, 9), 1, (4, 4, 4), (2, 2, 2))
    s = flex.sum(p)
    p = p / s

    # Copy profile
    c = add_poisson_noise(100 * p)
    b = flex.double(flex.grid(9, 9, 9), 0)
    m = flex.bool(flex.grid(9,9,9), True)

    # Fit
    fit = ProfileFitter(c, b, m, p)
    I = fit.intensity()
    V = fit.variance()
    assert fit.niter() < fit.maxiter()

    # Test intensity is the same
    eps = 1e-3
    assert(abs(I[0] - flex.sum(c)) < eps)
    assert(abs(V[0] - I[0]) < eps)

    print 'OK'

  def tst_with_flat_background(self):

    from dials.algorithms.integration.fit import ProfileFitter
    from scitbx.array_family import flex
    from numpy.random import seed
    seed(0)

    # Create profile
    p = gaussian((9, 9, 9), 1, (4, 4, 4), (2, 2, 2))
    s = flex.sum(p)
    p = p / s

    # Copy profile
    c0 = add_poisson_noise(100 * p)
    b = flex.double(flex.grid(9, 9, 9), 10)
    m = flex.bool(flex.grid(9,9,9), True)
    b0 = add_poisson_noise(b)
    c = c0 + b0

    # Fit
    fit = ProfileFitter(c, b, m, p)
    I = fit.intensity()
    V = fit.variance()
    assert fit.niter() < fit.maxiter()

    Iknown = 201.67417836585147
    Vknown = 7491.6743173001205

    # Test intensity is the same
    eps = 1e-3
    assert(abs(I[0] - Iknown) < eps)
    assert(abs(V[0] - Vknown) < eps)

    print 'OK'

  def tst_identical_partial(self):

    from dials.algorithms.integration.fit import ProfileFitter
    from scitbx.array_family import flex
    from numpy.random import seed
    seed(0)

    # Create profile
    p = gaussian((9, 9, 9), 1, (4, 4, 4), (2, 2, 2))
    s = flex.sum(p)
    p = p / s

    # Copy profile
    c = p.deep_copy()
    b = flex.double(flex.grid(9, 9, 9), 0)
    m = flex.bool(flex.grid(9,9,9), True)

    # Get the partial profiles
    pp = p[0:5,:,:]
    mp = m[0:5,:,:]
    cp = c[0:5,:,:]
    bp = b[0:5,:,:]

    # Fit
    fit = ProfileFitter(cp, bp, mp, pp)
    I = fit.intensity()
    V = fit.variance()
    assert fit.niter() < fit.maxiter()

    # Test intensity is the same
    eps = 1e-7
    assert(abs(I[0] - flex.sum(p)) < eps)
    assert(abs(V[0] - flex.sum(p)) < eps)

    print 'OK'

  def tst_with_no_background_partial(self):

    from dials.algorithms.integration.fit import ProfileFitter
    from scitbx.array_family import flex
    from numpy.random import seed
    seed(0)

    # Create profile
    p = gaussian((9, 9, 9), 1, (4, 4, 4), (2, 2, 2))
    s = flex.sum(p)
    p = p / s

    # Copy profile
    c = add_poisson_noise(100 * p)
    b = flex.double(flex.grid(9, 9, 9), 0)
    m = flex.bool(flex.grid(9,9,9), True)

    # Get the partial profiles
    pp = p[0:5,:,:]
    mp = m[0:5,:,:]
    cp = c[0:5,:,:]
    bp = b[0:5,:,:]

    # Fit
    fit = ProfileFitter(cp, bp, mp, pp)
    I = fit.intensity()
    V = fit.variance()
    assert fit.niter() < fit.maxiter()

    Iknown = 94.67151440319306
    Vknown = 94.67151440319306

    # Test intensity is the same
    eps = 1e-7
    assert(abs(I[0] - Iknown) < eps)
    assert(abs(V[0] - Vknown) < eps)

    print 'OK'

  def tst_with_flat_background_partial(self):

    from dials.algorithms.integration.fit import ProfileFitter
    from scitbx.array_family import flex
    from numpy.random import seed
    seed(0)

    # Create profile
    p = gaussian((9, 9, 9), 1, (4, 4, 4), (2, 2, 2))
    s = flex.sum(p)
    p = p / s

    # Copy profile
    c0 = add_poisson_noise(100 * p)
    b = flex.double(flex.grid(9, 9, 9), 1)
    m = flex.bool(flex.grid(9,9,9), True)
    c = c0 + add_poisson_noise(b)

    # Get the partial profiles
    pp = p[0:5,:,:]
    mp = m[0:5,:,:]
    cp = c[0:5,:,:]
    bp = b[0:5,:,:]
    c0p = c0[0:5,:,:]

    # Fit
    fit = ProfileFitter(cp, bp, mp, pp)
    I = fit.intensity()
    V = fit.variance()
    assert fit.niter() < fit.maxiter()

    Iknown = 99.06932141277105
    Vknown = 504.06932141277105

    # Test intensity is the same
    eps = 1e-7
    assert(abs(I[0] - Iknown) < eps)
    assert(abs(V[0] - Vknown) < eps)

    print 'OK'

  def generate_3_profiles(self):
    from dials.array_family import flex
    p1 = gaussian((40, 9, 9), 1, (10.5, 4, 4), (2, 2, 2))
    p2 = gaussian((40, 9, 9), 1, (20.5, 4, 4), (2, 2, 2))
    p3 = gaussian((40, 9, 9), 1, (30.5, 4, 4), (2, 2, 2))
    p1 = p1 / flex.sum(p1)
    p2 = p2 / flex.sum(p2)
    p3 = p3 / flex.sum(p3)
    p1.reshape(flex.grid(1, 40, 9, 9))
    p2.reshape(flex.grid(1, 40, 9, 9))
    p3.reshape(flex.grid(1, 40, 9, 9))
    p = flex.double(flex.grid(3, 40, 9, 9))
    p[0:1,:,:,:] = p1
    p[1:2,:,:,:] = p2
    p[2:3,:,:,:] = p3
    return p

  def generate_7_profiles(self):
    from dials.array_family import flex
    p1 = gaussian((40, 40, 40), 1, (10.5, 20.5, 20.5), (2, 2, 2))
    p2 = gaussian((40, 40, 40), 1, (20.5, 20.5, 20.5), (2, 2, 2))
    p3 = gaussian((40, 40, 40), 1, (30.5, 20.5, 20.5), (2, 2, 2))
    p4 = gaussian((40, 40, 40), 1, (20.5, 10.5, 20.5), (2, 2, 2))
    p5 = gaussian((40, 40, 40), 1, (20.5, 30.5, 20.5), (2, 2, 2))
    p6 = gaussian((40, 40, 40), 1, (20.5, 20.5, 10.5), (2, 2, 2))
    p7 = gaussian((40, 40, 40), 1, (20.5, 20.5, 30.5), (2, 2, 2))
    p1 = p1 / flex.sum(p1)
    p2 = p2 / flex.sum(p2)
    p3 = p3 / flex.sum(p3)
    p4 = p4 / flex.sum(p4)
    p5 = p5 / flex.sum(p5)
    p6 = p6 / flex.sum(p6)
    p7 = p7 / flex.sum(p7)
    p1.reshape(flex.grid(1, 40, 40, 40))
    p2.reshape(flex.grid(1, 40, 40, 40))
    p3.reshape(flex.grid(1, 40, 40, 40))
    p4.reshape(flex.grid(1, 40, 40, 40))
    p5.reshape(flex.grid(1, 40, 40, 40))
    p6.reshape(flex.grid(1, 40, 40, 40))
    p7.reshape(flex.grid(1, 40, 40, 40))
    p = flex.double(flex.grid(7, 40, 40, 40))
    p[0:1,:,:,:] = p1
    p[1:2,:,:,:] = p2
    p[2:3,:,:,:] = p3
    p[3:4,:,:,:] = p4
    p[4:5,:,:,:] = p5
    p[5:6,:,:,:] = p6
    p[6:7,:,:,:] = p7
    return p

  def tst_deconvolve_zero(self):
    from dials.algorithms.integration.fit import ProfileFitter
    from scitbx.array_family import flex
    from numpy.random import seed
    seed(0)

    I0 = [1000, 2000, 3000]

    # Create profile
    p = self.generate_3_profiles()


    # Copy profile
    c = flex.double(flex.grid(40, 9, 9), 0)
    b = flex.double(flex.grid(40, 9, 9), 0)
    m = flex.bool(flex.grid(40,9,9), True)

    # Fit
    passed = False
    try:
      fit = ProfileFitter(c, b, m, p)
      I = fit.intensity()
      V = fit.variance()
      passed = True
    except Exception:
      pass
    if passed:
      assert False, "This should fail"

    print 'OK'

  def tst_deconvolve_3_with_no_background(self):
    from dials.algorithms.integration.fit import ProfileFitter
    from scitbx.array_family import flex
    from numpy.random import seed
    seed(0)

    I0 = [1000, 2000, 3000]

    # Create profile
    p = self.generate_3_profiles()

    Ical = [[],[],[]]
    for it in range(1):

      # Copy profile
      c = flex.double(flex.grid(40, 9, 9))
      for i in range(p.all()[0]):
        pp = p[i:i+1,:,:,:]
        pp.reshape(c.accessor())
        cc = add_poisson_noise(I0[i] * pp)
        c += cc
      b = flex.double(flex.grid(40, 9, 9), 0)
      m = flex.bool(flex.grid(40,9,9), True)

      # Fit
      fit = ProfileFitter(c, b, m, p)
      I = fit.intensity()
      V = fit.variance()
      assert fit.niter() < fit.maxiter()

      for i in range(3):
        Ical[i].append(I[i])

    for i in range(3):
      Ical[i] = sum(Ical[i]) / len(Ical[i])

    Iknown = [1048.3221116842406, 1920.9035376774107, 2938.7743506383745]

    # Test intensity is the same
    eps = 1e-7
    for i in range(3):
      assert(abs(I[i] - Iknown[i]) < eps)
      assert(abs(V[i] - Iknown[i]) < eps)

    print 'OK'

  def tst_deconvolve_3_with_flat_background(self):
    from dials.algorithms.integration.fit import ProfileFitter
    from scitbx.array_family import flex
    from numpy.random import seed
    seed(0)

    I0 = [1000, 2000, 3000]

    # Create profile
    p = self.generate_3_profiles()

    Ical = [[],[],[]]
    for it in range(1):

      # Copy profile
      c = flex.double(flex.grid(40, 9, 9))
      for i in range(p.all()[0]):
        pp = p[i:i+1,:,:,:]
        pp.reshape(c.accessor())
        cc = add_poisson_noise(I0[i] * pp)
        c += cc
      b = flex.double(flex.grid(40, 9, 9), 1)
      m = flex.bool(flex.grid(40,9,9), True)
      c += add_poisson_noise(b)

      # Fit
      fit = ProfileFitter(c, b, m, p)
      I = fit.intensity()
      V = fit.variance()
      assert fit.niter() < fit.maxiter()

      for i in range(3):
        Ical[i].append(I[i])

    for i in range(3):
      Ical[i] = sum(Ical[i]) / len(Ical[i])

    Iknown = [1030.7018033805357, 1948.7152700952695, 2972.983204218213]
    Vknown = [4270.701803380534, 5188.715270095279, 6212.983204218214]

    # Test intensity is the same
    eps = 1e-7
    for i in range(3):
      assert(abs(I[i] - Iknown[i]) < eps)
      assert(abs(V[i] - Vknown[i]) < eps)

    print 'OK'

  def tst_deconvolve_7_with_no_background(self):
    from dials.algorithms.integration.fit import ProfileFitter
    from scitbx.array_family import flex
    from numpy.random import seed
    seed(0)

    I0 = [1000, 1500, 2000, 2500, 3000, 3500, 4000]

    # Create profile
    p = self.generate_7_profiles()

    Ical = [[],[],[],[],[],[],[]]
    for it in range(1):

      # Copy profile
      c = flex.double(flex.grid(40, 40, 40))
      for i in range(p.all()[0]):
        pp = p[i:i+1,:,:,:]
        pp.reshape(c.accessor())
        cc = add_poisson_noise(I0[i] * pp)
        c += cc
      b = flex.double(flex.grid(40, 40, 40), 0)
      m = flex.bool(flex.grid(40,40,40), True)

      # Fit
      fit = ProfileFitter(c, b, m, p)
      I = fit.intensity()
      V = fit.variance()
      assert fit.niter() < fit.maxiter()

      for i in range(7):
        Ical[i].append(I[i])

    for i in range(7):
      Ical[i] = sum(Ical[i]) / len(Ical[i])

    Iknown = [1012.4193633595916, 1472.3322059495797, 2072.136413825826,
              2486.4902438469253, 3012.6132119521126, 3409.2053517072773,
              3952.803209358826]

    # Test intensity is the same
    eps = 1e-7
    for i in range(7):
      assert(abs(I[i] - Iknown[i]) < eps)
      assert(abs(V[i] - Iknown[i]) < eps)

    print 'OK'

  def tst_deconvolve_7_with_flat_background(self):
    from dials.algorithms.integration.fit import ProfileFitter
    from scitbx.array_family import flex
    from numpy.random import seed
    seed(0)

    I0 = [1000, 1500, 2000, 2500, 3000, 3500, 4000]

    # Create profile
    p = self.generate_7_profiles()

    Ical = [[],[],[],[],[],[],[]]
    for it in range(1):

      # Copy profile
      c = flex.double(flex.grid(40, 40, 40))
      for i in range(p.all()[0]):
        pp = p[i:i+1,:,:,:]
        pp.reshape(c.accessor())
        cc = add_poisson_noise(I0[i] * pp)
        c += cc
      b = flex.double(flex.grid(40, 40, 40), 1)
      m = flex.bool(flex.grid(40,40,40), True)
      c += add_poisson_noise(b)

      # Fit
      fit = ProfileFitter(c, b, m, p)
      I = fit.intensity()
      V = fit.variance()
      assert fit.niter() < fit.maxiter()

      for i in range(7):
        Ical[i].append(I[i])

    for i in range(7):
      Ical[i] = sum(Ical[i]) / len(Ical[i])

    Iknown = [1042.4904898451261, 1447.6413897568257, 2053.2524102387893,
              2485.2789614429703, 3063.867598583333, 3445.9230910807582,
              3977.6744651425543]
    Vknown = [65042.49048984377, 65447.6413897533, 66053.25241023625,
              66485.27896144328, 67063.86759858252, 67445.92309107889,
              67977.67446514373]

    # Test intensity is the same
    eps = 1e-7
    for i in range(7):
      assert(abs(I[i] - Iknown[i]) < eps)
      assert(abs(V[i] - Vknown[i]) < eps)

    print 'OK'


if __name__ == '__main__':

  from dials.test import cd_auto
  with cd_auto(__file__):
    test = Test()
    test.run()

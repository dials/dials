from __future__ import absolute_import, division
from dials.test.algorithms.integration.profile.tst_profile_helpers import gaussian

class Test(object):

  def __init__(self):
    pass

  def run(self):
    self.tst_identical()
    self.tst_scaled()
    self.tst_with_flat_background()
    self.tst_with_noisy_flat_background()
    self.tst_identical_partial()
    self.tst_scaled_partial()
    self.tst_with_flat_background_partial()
    self.tst_with_noisy_flat_background_partial()

  def tst_identical(self):

    from dials.algorithms.integration.fit import fit_profile
    from scitbx.array_family import flex

    # Create profile
    p = gaussian((9, 9, 9), 1, (4, 4, 4), (2, 2, 2))
    s = flex.sum(p)
    p = p / s

    # Copy profile
    c = p.deep_copy()
    b = flex.double(flex.grid(9, 9, 9), 0)
    m = flex.bool(flex.grid(9,9,9), True)

    # Fit
    fit = fit_profile(p, m, c, b)
    I = fit.intensity()
    V = fit.variance()

    # Test intensity is the same
    eps = 1e-7
    assert(abs(I - flex.sum(p)) < eps)
    assert(abs(V - I) < eps)

    print 'OK'

  def tst_scaled(self):

    from dials.algorithms.integration.fit import fit_profile
    from scitbx.array_family import flex

    # Create profile
    p = gaussian((9, 9, 9), 1, (4, 4, 4), (2, 2, 2))
    s = flex.sum(p)
    p = p / s

    # Copy profile
    c = gaussian((9, 9, 9), 1, (4, 4, 4), (2, 2, 2))
    b = flex.double(flex.grid(9, 9, 9), 0)
    m = flex.bool(flex.grid(9,9,9), True)

    # Fit
    fit = fit_profile(p, m, c, b)
    I = fit.intensity()
    V = fit.variance()

    # Test intensity is the same
    eps = 1e-7
    assert(abs(I - flex.sum(c)) < eps)
    assert(abs(V - I) < eps)

    print 'OK'

  def tst_with_flat_background(self):

    from dials.algorithms.integration.fit import fit_profile
    from scitbx.array_family import flex

    # Create profile
    p = gaussian((9, 9, 9), 1, (4, 4, 4), (2, 2, 2))
    s = flex.sum(p)
    p = p / s

    # Copy profile
    c0 = gaussian((9, 9, 9), 1, (4, 4, 4), (2, 2, 2))
    b = flex.double(flex.grid(9, 9, 9), 5)
    m = flex.bool(flex.grid(9,9,9), True)
    c = c0 + b

    # Fit
    fit = fit_profile(p, m, c, b)
    I = fit.intensity()
    V = fit.variance()

    # Test intensity is the same
    eps = 1e-7
    assert(abs(I - flex.sum(c0)) < eps)
    assert(abs(V - (flex.sum(c0) + flex.sum(b))) < eps)

    print 'OK'

  def tst_with_noisy_flat_background(self):

    from dials.algorithms.integration.fit import fit_profile
    from scitbx.array_family import flex

    # Create profile
    p = gaussian((9, 9, 9), 1, (4, 4, 4), (2, 2, 2))
    s = flex.sum(p)
    p = p / s

    # Copy profile
    c0 = gaussian((9, 9, 9), 1, (4, 4, 4), (2, 2, 2))
    n = flex.random_double(9 * 9 * 9)
    b = flex.double(flex.grid(9, 9, 9), 5) + n
    m = flex.bool(flex.grid(9,9,9), True)
    c = c0 + b

    # Fit
    fit = fit_profile(p, m, c, b)
    I = fit.intensity()
    V = fit.variance()

    # Test intensity is the same
    eps = 1e-7
    assert(abs(I - flex.sum(c0)) < eps)
    assert(abs(V - (flex.sum(c0) + flex.sum(b))) < eps)

    print 'OK'

  def tst_identical_partial(self):

    from dials.algorithms.integration.fit import fit_profile
    from scitbx.array_family import flex

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
    fit = fit_profile(pp, mp, cp, bp)
    I = fit.intensity()
    V = fit.variance()

    # Test intensity is the same
    eps = 1e-7
    assert(abs(I - flex.sum(p)) < eps)
    assert(abs(V - flex.sum(pp)) < eps)

    print 'OK'

  def tst_scaled_partial(self):

    from dials.algorithms.integration.fit import fit_profile
    from scitbx.array_family import flex

    # Create profile
    p = gaussian((9, 9, 9), 1, (4, 4, 4), (2, 2, 2))
    s = flex.sum(p)
    p = p / s

    # Copy profile
    c = gaussian((9, 9, 9), 1, (4, 4, 4), (2, 2, 2))
    b = flex.double(flex.grid(9, 9, 9), 0)
    m = flex.bool(flex.grid(9,9,9), True)

    # Get the partial profiles
    pp = p[0:5,:,:]
    mp = m[0:5,:,:]
    cp = c[0:5,:,:]
    bp = b[0:5,:,:]

    # Fit
    fit = fit_profile(pp, mp, cp, bp)
    I = fit.intensity()
    V = fit.variance()

    # Test intensity is the same
    eps = 1e-7
    assert(abs(I - flex.sum(c)) < eps)
    assert(abs(V - flex.sum(cp)) < eps)

    print 'OK'

  def tst_with_flat_background_partial(self):

    from dials.algorithms.integration.fit import fit_profile
    from scitbx.array_family import flex

    # Create profile
    p = gaussian((9, 9, 9), 1, (4, 4, 4), (2, 2, 2))
    s = flex.sum(p)
    p = p / s

    # Copy profile
    c0 = gaussian((9, 9, 9), 1, (4, 4, 4), (2, 2, 2))
    b = flex.double(flex.grid(9, 9, 9), 5)
    m = flex.bool(flex.grid(9,9,9), True)
    c = c0 + b

    # Get the partial profiles
    pp = p[0:5,:,:]
    mp = m[0:5,:,:]
    cp = c[0:5,:,:]
    bp = b[0:5,:,:]
    c0p = c0[0:5,:,:]

    # Fit
    fit = fit_profile(pp, mp, cp, bp)
    I = fit.intensity()
    V = fit.variance()

    # Test intensity is the same
    eps = 1e-7
    assert(abs(I - flex.sum(c0)) < eps)
    assert(abs(V - (flex.sum(c0p) + flex.sum(bp))) < eps)

    print 'OK'

  def tst_with_noisy_flat_background_partial(self):

    from dials.algorithms.integration.fit import fit_profile
    from scitbx.array_family import flex

    # Create profile
    p = gaussian((9, 9, 9), 1, (4, 4, 4), (2, 2, 2))
    s = flex.sum(p)
    p = p / s

    # Copy profile
    c0 = gaussian((9, 9, 9), 1, (4, 4, 4), (2, 2, 2))
    n = flex.random_double(9 * 9 * 9)
    b = flex.double(flex.grid(9, 9, 9), 5) + n
    m = flex.bool(flex.grid(9,9,9), True)
    c = c0 + b

    # Get the partial profiles
    pp = p[0:5,:,:]
    mp = m[0:5,:,:]
    cp = c[0:5,:,:]
    bp = b[0:5,:,:]
    c0p = c0[0:5,:,:]

    # Fit
    fit = fit_profile(pp, mp, cp, bp)
    I = fit.intensity()
    V = fit.variance()

    # Test intensity is the same
    eps = 1e-7
    assert(abs(I - flex.sum(c0)) < eps)
    assert(abs(V - (flex.sum(c0p) + flex.sum(bp))) < eps)

    print 'OK'

if __name__ == '__main__':

  from dials.test import cd_auto
  with cd_auto(__file__):
    test = Test()
    test.run()

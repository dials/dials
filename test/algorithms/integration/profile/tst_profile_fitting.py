from __future__ import division

class Test(object):

    def __init__(self):
        pass

    def run(self):
        self.tst_identical()
        self.tst_scaled()

    def tst_identical(self):

        from dials.algorithms.integration.profile import ProfileFitting
        from scitbx.array_family import flex
        from tst_profile_helpers import gaussian

        # Create profile
        p = gaussian((9, 9, 9), 1, (4, 4, 4), (2, 2, 2))
        s = flex.sum(p)
        p = p / s

        # Copy profile
        c = p.deep_copy()
        b = flex.double(flex.grid(9, 9, 9), 0)

        # Fit
        fit = ProfileFitting(p, c, b)
        I = fit.intensity()
        V = fit.variance()

        # Test intensity is the same
        eps = 1e-7
        assert(abs(I - flex.sum(p)) < eps)
        assert(abs(V - I) < eps)

        print 'OK'

    def tst_scaled(self):

        from dials.algorithms.integration.profile import ProfileFitting
        from scitbx.array_family import flex
        from tst_profile_helpers import gaussian

        # Create profile
        p = gaussian((9, 9, 9), 1, (4, 4, 4), (2, 2, 2))
        s = flex.sum(p)
        p = p / s

        # Copy profile
        c = gaussian((9, 9, 9), 1, (4, 4, 4), (2, 2, 2))
        b = flex.double(flex.grid(9, 9, 9), 0)

        # Fit
        fit = ProfileFitting(p, c, b)
        I = fit.intensity()
        V = fit.variance()

        # Test intensity is the same
        eps = 1e-7
        assert(abs(I - flex.sum(c)) < eps)
        assert(abs(V - I) < eps)

        print 'OK'

if __name__ == '__main__':
    test = Test()
    test.run()

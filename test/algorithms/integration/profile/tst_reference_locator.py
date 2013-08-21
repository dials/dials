from __future__ import division

class Test(object):

    def __init__(self):
        from scitbx.array_family import flex
        from dials.algorithms.integration.profile import XdsCircleSampler
        self.sampler = XdsCircleSampler((1000, 1000, 10), 2)

        self.profiles = flex.double(flex.grid(
          len(self.sampler),
          9, 9, 9))

        for i in range(len(self.sampler) * 9 * 9 * 9):
            self.profiles[i] = i

    def run(self):

        from dials.algorithms.integration.profile import \
            XdsCircleReferenceLocator
        import tempfile
        import cPickle as pickle

        # Create the reference locator
        locator = XdsCircleReferenceLocator(self.profiles, self.sampler)

        # Test basic access
        assert(locator.size() == len(self.sampler))
        profiles = locator.profile()
        assert(profiles.all() == self.profiles.all())

        print 'OK'

        # Test access by index
        for i in range(locator.size()):
            profile = locator.profile(i)
            assert(profile.all() == (9, 9, 9))

            # Check the values
            for j in range(len(profile)):
                assert(profile[j] == j + i*9*9*9)

        print 'OK'

        # Dump the locator to pickle
        tf = tempfile.TemporaryFile()
        pickle.dump(locator, tf)
        tf.flush()
        tf.seek(0)
        locator2 = pickle.load(tf)

        # Check the profiles are the same
        p1 = locator.profile()
        p2 = locator2.profile()
        assert(len(p1) == len(p2))
        assert(p1.all() == p2.all())
        eps = 1e-7
        for i in range(len(p1)):
            assert((p1[i] - p2[i]) <= eps)

        # Check the sampler is the same
        s1 = locator.sampler()
        s2 = locator2.sampler()
        assert(s1.volume_size() == s2.volume_size())
        assert(s1.num_z() == s2.num_z())

        # test passed
        print 'OK'

if __name__ == '__main__':
    test = Test()
    test.run()

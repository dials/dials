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


if __name__ == '__main__':
    test = Test()
    test.run()

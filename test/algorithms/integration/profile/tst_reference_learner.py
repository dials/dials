
class Test(object):

    def __init__(self):

        from dials.algorithms.integration.profile import XdsCircleSampler

        width = 1000
        height = 1000
        depth = 10
        nz = 2
        self.sampler = XdsCircleSampler((width, height, depth), nz)

        self.grid_size = (9, 9, 9)
        self.threshold = 0.2

    def run(self):

        self.tst_with_identical_non_negative_profiles()
        self.tst_with_systematically_offset_profiles()

    def tst_with_identical_non_negative_profiles(self):

        from dials.algorithms.integration.profile import XdsCircleReferenceLearner

        # Generate identical non-negative profiles
        reflections, profile = self.generate_identical_non_negative_profiles()

        # Create the reference learner
        learner = XdsCircleReferenceLearner(self.sampler,
            self.grid_size, self.threshold)

        # Learn from the reflections
        learner.learn(reflections)

        # Get the reference locator
        locate = learner.locate()

        # Normalize the profile
        profile = self.normalize_profile(profile)

        # Check that all the reference profiles are the same
        eps = 1e-10
        for index in range(locate.size()):
            reference = locate.profile(index)
            for k in range(self.grid_size[2]):
                for j in range(self.grid_size[1]):
                    for i in range(self.grid_size[0]):
                        assert(abs(reference[k,j,i] - profile[k,j,i]) <= eps)

        print 'OK'

    def tst_with_systematically_offset_profiles(self):
        from tst_profile_helpers import gaussian
        from dials.algorithms.integration.profile import XdsCircleReferenceLearner
        from dials.algorithms.image.centroid import CentroidImage3d
        from scitbx import matrix
        from math import sqrt, pi, cos

        # Generate identical non-negative profiles
        reflections = self.generate_systematically_offset_profiles()

        # Create the reference learner
        learner = XdsCircleReferenceLearner(self.sampler,
            self.grid_size, self.threshold)

        # Learn from the reflections
        learner.learn(reflections)

        # Get the reference locator
        locate = learner.locate()

        # Check that all the reference profiles are the same
        eps = 1e-3
        x1 = []
        x2 = []
        for index in range(locate.size()):
            reference = locate.profile(index)
            coord = locate.coord(index)
            centroid = CentroidImage3d(reference)
            x1.append(coord[0])
            x2.append(centroid.mean()[0])

        Y = matrix.col((
            sum(x2),
            sum([x * y for x, y in zip(x1, x2)])))
        X = matrix.sqr((
            len(x1),
            sum(x1),
            sum(x1),
            sum([x * x for x in x1])))

        b = X.inverse() * Y

        print 'OK'

    def normalize_profile(self, profile):
        from scitbx.array_family import flex
        max_profile = flex.max(profile)
        threshold = self.threshold * max_profile
        sum_profile = 0.0
        for i in range(len(profile)):
            if profile[i] > threshold:
                sum_profile += profile[i]

        result = flex.double(flex.grid(profile.all()))
        for i in range(len(profile)):
            result[i] = profile[i] / sum_profile

        return result


    def generate_identical_non_negative_profiles(self):
        from dials.model.data import ReflectionList, Reflection
        from scitbx.array_family import flex
        from random import uniform
        from tst_profile_helpers import gaussian
        rlist = ReflectionList(1000)

        profile = gaussian(self.grid_size, 1000, (4, 4, 4), (1.5, 1.5, 1.5))

        for i in range(1000):
            x = uniform(0, 1000)
            y = uniform(0, 1000)
            z = uniform(0, 10)
            rlist[i].status = 0
            rlist[i].image_coord_px = (x, y)
            rlist[i].frame_number = z
            rlist[i].transformed_shoebox = profile.deep_copy()

        return rlist, profile

    def generate_systematically_offset_profiles(self):
        from dials.model.data import ReflectionList, Reflection
        from scitbx.array_family import flex
        from random import uniform
        from tst_profile_helpers import gaussian
        rlist = ReflectionList(1000)

        for i in range(1000):
            x = uniform(0, 1000)
            y = uniform(0, 1000)
            z = uniform(0, 10)

            offset = -4.5  + 9 * x / 1000.0

            profile = gaussian(self.grid_size, 1000,
                (4 + offset, 4, 4), (1.5, 1.5, 1.5))
            rlist[i].status = 0
            rlist[i].image_coord_px = (x, y)
            rlist[i].frame_number = z
            rlist[i].transformed_shoebox = profile

        return rlist

if __name__ == '__main__':
    test = Test()
    test.run()

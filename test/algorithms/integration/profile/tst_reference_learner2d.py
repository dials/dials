from __future__ import division

class Test(object):

  def __init__(self):

    from dials.algorithms.integration.profile import GridSampler2D
    from dials.array_family import flex

    width = 1000
    height = 1000
    self.sampler = GridSampler2D((width, height), (5, 5))
    self.min_average = 10
    self.nrefl = 1000
    self.xsize = flex.size_t([9] * 5*5)
    self.ysize = flex.size_t([9] * 5*5)

  def run(self):

    self.tst_with_identical_non_negative_profiles()
    self.tst_with_systematically_offset_profiles()
    # self.tst_with_single_profile()

  def tst_with_identical_non_negative_profiles(self):

    from dials.algorithms.integration.profile import ReferenceLearner2D
    from scitbx.array_family import flex
    from dials.algorithms.shoebox import MaskCode

    # Generate identical non-negative profiles
    xyz, profiles, profile = self.generate_identical_non_negative_profiles()

    # Create the reference learner
    learner = ReferenceLearner2D(self.sampler, self.xsize, self.ysize,
                                 self.min_average)

    # Learn from the reflections
    for i in range(1000):
      data = profiles[i]
      background = flex.double(flex.grid(data.all()), 0)
      mask = flex.int(flex.grid(data.all()),
                      MaskCode.Valid | MaskCode.Foreground)
      coord = xyz[i]
      learner.add(data, background, mask, coord)
    learner.finish()

    # Normalize the profile
    profile = self.normalize_profile(profile)

    # Check that all the reference profiles are the same
    eps = 1e-10
    for index in range(len(learner)):
      reference = learner.profile(index)
      for j in range(9):
        for i in range(9):
          assert(abs(reference[j,i] - profile[j,i]) <= eps)
      assert(abs(flex.sum(reference) - 1.0) <= eps)

    print 'OK'

  def tst_with_systematically_offset_profiles(self):
    from tst_profile_helpers import gaussian
    from dials.algorithms.integration.profile import ReferenceLearner2D
    from dials.algorithms.image.centroid import centroid_image
    from scitbx import matrix
    from math import sqrt, pi, cos
    from scitbx.array_family import flex
    from dials.algorithms.shoebox import MaskCode

    # Generate identical non-negative profiles
    xyz, profiles = self.generate_systematically_offset_profiles()

    # Create the reference learner
    learner = ReferenceLearner2D(self.sampler, self.xsize, self.ysize,
                                 self.min_average)

    # Learn from the reflections
    for i in range(1000):
      data = profiles[i]
      background = flex.double(flex.grid(data.all()), 0)
      mask = flex.int(flex.grid(data.all()),
                      MaskCode.Valid | MaskCode.Foreground)
      coord = xyz[i]
      learner.add(data, background, mask, coord)
    learner.finish()

    # Check that all the reference profiles are the same
    eps = 1e-3
    x1 = []
    x2 = []
    for index in range(len(learner)):
      reference = learner.profile(index)
      coord = learner.coord(index)
      centroid = centroid_image(reference)
      x1.append(coord[0])
      x2.append(centroid.mean()[0])

      assert(abs(flex.sum(reference) - 1.0) < 1e-7)

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

  # def tst_with_single_profile(self):

  #   from dials.algorithms.integration.profile import ReferenceLearner2D
  #   from scitbx.array_family import flex
  #   from dials.algorithms.shoebox import MaskCode

  #   # Generate identical non-negative profiles
  #   xyz, profiles, profile = self.generate_single_central_non_negative_profiles()

  #   # Create the reference learner
  #   learner = ReferenceLearner2D(self.sampler, self.xsize, self.ysize,
  #                                self.min_average)

  #   # Learn from the reflections
  #   for i in range(1):
  #     data = profiles[i]
  #     background = flex.double(flex.grid(data.all()), 0)
  #     mask = flex.int(flex.grid(data.all()),
  #                     MaskCode.Valid | MaskCode.Foreground)
  #     coord = xyz[i]
  #     learner.add(data, background, mask, coord)
  #   learner.finish()

  #   # Normalize the profile
  #   profile = self.normalize_profile(profile)
  #   zero = flex.double(profile.accessor(), 0)

  #   assert(len(reflections) == 1)
  #   coord = reflections[0]['xyzcal.px']
  #   ind = locate.indices(coord)

  #   nind = set(range(locate.size())).difference(ind)
  #   assert(len(nind) == 0)

  #   for i in ind:
  #     assert(locate.profile(i).as_1d().all_approx_equal(profile.as_1d()))

  #   print 'OK'


  def normalize_profile(self, profile):
    from scitbx.array_family import flex
    max_profile = flex.max(profile)
    # threshold = self.threshold * max_profile
    threshold = 0
    sum_profile = 0.0
    for i in range(len(profile)):
      if profile[i] > threshold:
        sum_profile += profile[i]
      else:
        profile[i] = 0.0

    result = flex.double(flex.grid(profile.all()))
    for i in range(len(profile)):
      result[i] = profile[i] / sum_profile

    return result


  def generate_single_central_non_negative_profiles(self):
    from dials.array_family import flex
    from random import uniform
    from tst_profile_helpers import gaussian

    profile = gaussian((9, 9), 1000, (4, 4), (1.5, 1.5))

    x = 500
    y = 500
    z = 5
    xyz = flex.vec3_double(1)
    xyz[0] = (x, y, z)
    profiles = [None]
    profiles[0] = profile.deep_copy()

    return xyz, profiles, profile


  def generate_identical_non_negative_profiles(self):
    from dials.array_family import flex
    from random import uniform
    from tst_profile_helpers import gaussian

    profile = gaussian((9, 9), 1000, (4, 4), (1.5, 1.5))

    xyz = flex.vec3_double(1000)
    profiles = [None for i in range(1000)]
    for i in range(1000):
      x = uniform(0, 1000)
      y = uniform(0, 1000)
      z = uniform(0, 10)
      xyz[i] = (x, y, z)
      profiles[i] = profile.deep_copy()
    return xyz, profiles, profile

  def generate_systematically_offset_profiles(self):
    from dials.array_family import flex
    from random import uniform
    from tst_profile_helpers import gaussian

    xyz = flex.vec3_double(1000)
    profiles = [None for i in range(1000)]
    for i in range(1000):
      x = uniform(0, 1000)
      y = uniform(0, 1000)
      z = uniform(0, 10)

      offset = -4.5  + 9 * x / 1000.0

      profile = gaussian((9, 9), 1000,
          (4 + offset, 4), (1.5, 1.5))
      xyz[i] = (x, y, z)
      profiles[i] = profile

    return xyz, profiles

if __name__ == '__main__':
  from dials.test import cd_auto
  with cd_auto(__file__):
    test = Test()
    test.run()

from __future__ import absolute_import, division
from dials.algorithms.profile_model.modeller import EmpiricalProfileModeller
from dials.test.algorithms.integration.profile.tst_profile_helpers import gaussian

class Modeller(EmpiricalProfileModeller):

  def model(self, reflections, profiles):
    from dials.array_family import flex
    indices = flex.size_t(range(len(self)))
    weights = flex.double([1.0]*len(self))
    for profile in profiles:
      self.add(indices, weights, profile)

class Test(object):

  def __init__(self):

    width = 1000
    height = 1000
    depth = 10
    self.n = 9
    self.grid_size = (9, 9, 9)
    self.threshold = 0.0

  def run(self):

    self.tst_with_identical_non_negative_profiles()
    self.tst_with_systematically_offset_profiles()

  def tst_with_identical_non_negative_profiles(self):

    from scitbx.array_family import flex

    # Generate identical non-negative profiles
    reflections, profiles, profile = self.generate_identical_non_negative_profiles()

    # Create the reference learner
    modeller = Modeller(self.n, self.grid_size, self.threshold)

    # Do the modelling
    modeller.model(reflections, profiles)
    modeller.finalize()

    # Normalize the profile
    profile = self.normalize_profile(profile)

    # Check that all the reference profiles are the same
    eps = 1e-10
    for index in range(len(modeller)):
      reference = modeller.data(index)
      for k in range(self.grid_size[2]):
        for j in range(self.grid_size[1]):
          for i in range(self.grid_size[0]):
            assert(abs(reference[k,j,i] - profile[k,j,i]) <= eps)
      assert(abs(flex.sum(reference) - 1.0) <= eps)

    print 'OK'

  def tst_with_systematically_offset_profiles(self):
    from scitbx.array_family import flex

    # Generate identical non-negative profiles
    reflections, profiles = self.generate_systematically_offset_profiles()

    # Create the reference learner
    modeller = Modeller(self.n, self.grid_size, self.threshold)

    # Do the modelling
    modeller.model(reflections, profiles)
    modeller.finalize()

    # Check that all the reference profiles are the same
    eps = 1e-10
    profile = None
    for index in range(len(modeller)):
      reference = modeller.data(index)
      if profile is not None:
        for k in range(self.grid_size[2]):
          for j in range(self.grid_size[1]):
            for i in range(self.grid_size[0]):
              assert(abs(reference[k,j,i] - profile[k,j,i]) <= eps)
      else:
        profile = reference
      assert(abs(flex.sum(reference) - 1.0) <= eps)

    print 'OK'


  def normalize_profile(self, profile):
    from scitbx.array_family import flex
    max_profile = flex.max(profile)
    threshold = self.threshold * max_profile
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
    rlist = flex.reflection_table(1)

    profile = gaussian(self.grid_size, 1000, (4, 4, 4), (1.5, 1.5, 1.5))

    x = 500
    y = 500
    z = 5
    xyz = flex.vec3_double(1)
    xyz[0] = (x, y, z)
    profiles = [ profile.deep_copy() ]
    rlist['xyzcal.px'] = xyz

    return rlist, profiles, profile


  def generate_identical_non_negative_profiles(self):
    from dials.array_family import flex
    from random import uniform
    rlist = flex.reflection_table(1000)

    profile = gaussian(self.grid_size, 1000, (4, 4, 4), (1.5, 1.5, 1.5))

    xyz = flex.vec3_double(1000)
    profiles = []
    for i in range(1000):
      x = uniform(0, 1000)
      y = uniform(0, 1000)
      z = uniform(0, 10)
      xyz[i] = (x, y, z)
      profiles.append(profile.deep_copy())
    rlist['xyzcal.px'] = xyz

    return rlist, profiles, profile

  def generate_systematically_offset_profiles(self):
    from dials.array_family import flex
    from random import uniform
    rlist = flex.reflection_table(1000)

    xyz = flex.vec3_double(1000)
    profiles = []
    for i in range(1000):
      x = uniform(0, 1000)
      y = uniform(0, 1000)
      z = uniform(0, 10)

      offset = -4.5  + 9 * x / 1000.0

      profile = gaussian(self.grid_size, 1000,
          (4 + offset, 4, 4), (1.5, 1.5, 1.5))
      xyz[i] = (x, y, z)
      profiles.append(profile)

    rlist['xyzcal.px'] = xyz
    return rlist, profiles

if __name__ == '__main__':
  from dials.test import cd_auto
  with cd_auto(__file__):
    test = Test()
    test.run()


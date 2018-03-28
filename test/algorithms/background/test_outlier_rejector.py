from __future__ import absolute_import, division, print_function

def generate_shoebox(size, mean, nforeground, ninvalid):
  from dials.algorithms.simulation.generate_test_reflections \
    import random_background_plane2
  from dials.array_family import flex
  from random import sample
  from dials.algorithms.shoebox import MaskCode
  data = flex.double(flex.grid(size), 0.0)
  mask = flex.int(flex.grid(size), MaskCode.Valid | MaskCode.Background)
  random_background_plane2(data, mean, 0, 0, 0)
  for i in sample(range(len(data)), ninvalid):
    mask[i] &= ~MaskCode.Valid
  for i in sample(range(len(data)), nforeground):
    mask[i] |= MaskCode.Foreground
    mask[i] &= ~MaskCode.Background
  return data, mask

def assert_basic_mask_is_correct(mask, ninvalid, nforeground):
  from scitbx.array_family import flex
  from dials.algorithms.shoebox import MaskCode
  invalid = flex.size_t(
    i for i in range(len(mask)) if not mask[i] & MaskCode.Valid)
  foreground = flex.size_t(
    i for i in range(len(mask)) if mask[i] & MaskCode.Foreground)
  background = flex.size_t(
    i for i in range(len(mask)) if mask[i] & MaskCode.Background)
  background_used = flex.size_t(
    i for i in range(len(mask)) if mask[i] & MaskCode.BackgroundUsed)
  background_valid = flex.size_t(
    i for i in range(len(mask))
      if mask[i] & MaskCode.Valid and mask[i] & MaskCode.Background)
  assert(len(invalid) == ninvalid)
  assert(len(foreground) == nforeground)
  assert(len(background) == len(mask) - len(foreground))
  assert(len(background_used) < len(background))
  assert(len(set(background).intersection(foreground)) == 0)
  assert(len(set(background).intersection(background_used)) == len(background_used))
  return invalid, foreground, background, background_used, background_valid

def test_truncated():
  from dials.algorithms.background.simple import TruncatedOutlierRejector
  lower = 0.01
  upper = 0.01
  reject = TruncatedOutlierRejector(lower, upper)
  size = (9,9,9)
  ninvalid = 5
  nforeground = 20
  mean = 20

  def assert_is_correct(data, mask):
    from dials.array_family import flex
    invalid, foreground, background, background_used, background_valid  = \
      assert_basic_mask_is_correct(mask, ninvalid, nforeground)
    ntot = len(background_valid)
    i0 = int(lower * ntot / 2.0)
    i1 = int(upper * ntot / 2.0)
    nexp = ntot - i1 - i0
    assert len(background_used) == nexp
    index_exp = sorted(mask.select(background_valid), key=lambda x:data[x])[i0:i1]
    assert all(ii == jj for ii, jj in zip(index_exp, background_used))

  for i in range(10):
    data, mask = generate_shoebox(
        size,
        mean,
        nforeground,
        ninvalid)
    reject(data, mask)
    assert_is_correct(data, mask)

def test_nsigma():
  from dials.algorithms.background.simple import NSigmaOutlierRejector
  lower = 3.0
  upper = 3.0
  reject = NSigmaOutlierRejector(lower, upper)
  size = (9,9,9)
  ninvalid = 5
  nforeground = 20
  mean = 20

  def assert_is_correct(data, mask):
    from scitbx.array_family import flex
    invalid, foreground, background, background_used, background_valid = \
      assert_basic_mask_is_correct(mask, ninvalid, nforeground)
    subdata = data.select(background_valid)
    mv = flex.mean_and_variance(subdata)
    m = mv.mean()
    s = mv.unweighted_sample_standard_deviation()
    p0 = m - lower * s
    p1 = m + upper * s
    mask = (subdata >= p0) & (subdata <= p1)
    exp = background_valid.select(mask)
    assert len(exp) == len(background_used)
    assert all(ii == jj for ii, jj in zip(exp, background_used))

  for i in range(10):
    data, mask = generate_shoebox(
      size,
      mean,
      nforeground,
      ninvalid)
    reject(data, mask)
    assert_is_correct(data, mask)

def test_normal():
  from dials.algorithms.background.simple import NormalOutlierRejector
  min_data = 10
  reject = NormalOutlierRejector(min_data)
  size = (9,9,9)
  ninvalid = 5
  nforeground = 20
  mean = 20

  def assert_is_correct(data, mask):
    invalid, foreground, background, background_used, background_valid = \
      assert_basic_mask_is_correct(mask, ninvalid, nforeground)

  for i in range(10):
    data, mask = generate_shoebox(
      size,
      mean,
      nforeground,
      ninvalid)
    reject(data, mask)
    assert_is_correct(data, mask)

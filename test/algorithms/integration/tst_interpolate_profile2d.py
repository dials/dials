

from __future__ import absolute_import, division

class Test(object):

  def __init__(self):
    pass

  def run(self):
    from dials.algorithms.integration import interpolate_profile2d
    from dials.array_family import flex
    eps = 1e-7
    data = flex.double([
      [0, 0, 1, 0, 0],
      [0, 1, 2, 1, 0],
      [1, 2, 3, 2, 1],
      [0, 1, 2, 1, 0],
      [0, 0, 1, 0, 0]])

    x, y = (0.0, 0.0)
    result = interpolate_profile2d(x, y, data)
    exp = flex.double([
      [0, 0, 1, 0, 0],
      [0, 1, 2, 1, 0],
      [1, 2, 3, 2, 1],
      [0, 1, 2, 1, 0],
      [0, 0, 1, 0, 0]])
    for x1, x2 in zip(result, exp):
      assert(abs(x1 - x2) < eps)

    x, y = (-0.5, 0.0)
    result = interpolate_profile2d(x, y, data)
    exp = flex.double([
      [0,   0, 0.5, 0.5,   0],
      [0, 0.5, 1.5, 1.5, 0.5],
      [1, 1.5, 2.5, 2.5, 1.5],
      [0, 0.5, 1.5, 1.5, 0.5],
      [0,   0, 0.5, 0.5,   0]])
    for x1, x2 in zip(result, exp):
      assert(abs(x1 - x2) < eps)

    x, y = (0.0, -0.5)
    result = interpolate_profile2d(x, y, data)
    exp = flex.double([
      [  0,   0,   1,   0,   0],
      [  0, 0.5, 1.5, 0.5,   0],
      [0.5, 1.5, 2.5, 1.5, 0.5],
      [0.5, 1.5, 2.5, 1.5, 0.5],
      [  0, 0.5, 1.5, 0.5,   0]])
    for x1, x2 in zip(result, exp):
      assert(abs(x1 - x2) < eps)

    x, y = (0.5, 0.0)
    result = interpolate_profile2d(x, y, data)
    exp = flex.double([
      [  0, 0.5, 0.5,   0, 0],
      [0.5, 1.5, 1.5, 0.5, 0],
      [1.5, 2.5, 2.5, 1.5, 1],
      [0.5, 1.5, 1.5, 0.5, 0],
      [  0, 0.5, 0.5,   0, 0]])
    for x1, x2 in zip(result, exp):
      assert(abs(x1 - x2) < eps)

    x, y = (0.0, 0.5)
    result = interpolate_profile2d(x, y, data)
    exp = flex.double([
      [  0, 0.5, 1.5, 0.5,   0],
      [0.5, 1.5, 2.5, 1.5, 0.5],
      [0.5, 1.5, 2.5, 1.5, 0.5],
      [  0, 0.5, 1.5, 0.5,   0],
      [  0,   0,   1,   0,   0]])
    for x1, x2 in zip(result, exp):
      assert(abs(x1 - x2) < eps)

    x, y = -0.5, 0.5
    result = interpolate_profile2d(x, y, data)
    exp = flex.double([
      [  0, 0.25,   1,   1, 0.25],
      [0.5,    1,   2,   2,    1],
      [0.5,    1,   2,   2,    1],
      [  0, 0.25,   1,   1, 0.25],
      [  0,    0, 0.5, 0.5,    0]])
    for x1, x2 in zip(result, exp):
      assert(abs(x1 - x2) < eps)

    print 'OK'

if __name__ == '__main__':
  test = Test()
  test.run()

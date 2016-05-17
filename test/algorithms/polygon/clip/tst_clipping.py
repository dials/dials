from __future__ import division

def point_in_polygon(point, poly):
  inside = False
  j = len(poly) - 1
  for i in range(len(poly)):
    if (((poly[i][1] > point[1]) != (poly[j][1] > point[1])) and
        (point[0] < (poly[j][0] - poly[i][0]) * (point[1] - poly[i][1]) /
                    (poly[j][1] - poly[i][1]) + poly[i][0])):
      inside = not inside
  return inside

def generate_intersecting(subject_size=None, target_size=None):
  from random import randint

  if subject_size is None:
    subject_size = randint(3, 10)

  if target_size is None:
    target_size = 4

  bbox = (0.0, 10.0, 0.0, 10.0)
  subject = generate_polygon(subject_size, bbox)
  clip = generate_polygon(target_size, bbox)

  return subject, clip

def generate_non_intersecting(subject_size=None, target_size=None):
  from random import randint

  if subject_size is None:
    subject_size = randint(3, 10)

  if target_size is None:
    target_size = 4

  bbox = (0.0, 10.0, 0.0, 10.0)
  subject = generate_polygon(subject_size, bbox)
  target = generate_polygon(target_size, bbox)

  offset = [
      lambda: [(x + 10.0, y) for x, y in subject],
      lambda: [(x - 10.0, y) for x, y in subject],
      lambda: [(x, y + 10.0) for x, y in subject],
      lambda: [(x, y - 10.0) for x, y in subject],
      lambda: [(x + 10.0, y - 10.0) for x, y in subject],
      lambda: [(x - 10.0, y - 10.0) for x, y in subject],
      lambda: [(x + 10.0, y + 10.0) for x, y in subject],
      lambda: [(x - 10.0, y + 10.0) for x, y in subject],
  ]

  subject = offset[randint(0, 7)]()
  return subject, target

def generate_polygon(nvert, box):
  from random import uniform
  from math import pi, sin, cos
  xc = (box[0] + box[1]) / 2
  yc = (box[2] + box[3]) / 2
  maxr = min([xc, yc])
  minr = 0.1 * maxr
  v = []
  angle = 2 * pi
  dt = 2 * pi / nvert
  for i in range(nvert):
    r = uniform(minr, maxr)
    x = r * cos(angle)
    y = r * sin(angle)
    v.append((x, y))
    angle = angle + dt

  return v

class TestSimpleWithConvex(object):

  def __init__(self):
    pass

  def __call__(self):
    self.tst_intersecting()
    self.tst_non_intersecting()

  def tst_intersecting(self):
    from dials.algorithms.polygon import clip
    from scitbx.array_family import flex
    for i in range(10000):

      # Generate intersecting polygons
      subject, target = generate_intersecting()

      # Do the clipping
      result = clip.simple_with_convex(
          flex.vec2_double(subject),
          flex.vec2_double(target))

      # Ensure we have roughly valid number of vertices
      assert(len(result) >= 3)
      assert(len(result) >= min([len(subject), len(target)]))

#            for v in result:
#                assert(point_in_polygon(v, clip))

    print 'OK'

  def tst_non_intersecting(self):
    from dials.algorithms.polygon import clip
    from scitbx.array_family import flex

    for i in range(10000):

      # Generate nonintersecting polygons
      subject, target = generate_non_intersecting()

      # Do the clipping
      result = clip.simple_with_convex(
          flex.vec2_double(subject),
          flex.vec2_double(target))

      # Ensure we no vertices
      assert(len(result) == 0)

    print 'OK'


class TestSimpleWithRect(object):

  def __init__(self):
    pass

  def __call__(self):
    self.tst_intersecting()
    #self.tst_non_intersecting()

  def tst_intersecting(self):
    from dials.algorithms.polygon import clip
    from scitbx.array_family import flex
    for i in range(10000):

      # Generate intersecting polygons
      subject, target = generate_intersecting(target_size=2)
      rect = ((0, 0), (10, 10))

      # Do the clipping
      result = clip.simple_with_rect(
          flex.vec2_double(subject), rect)

      # Ensure we have roughly valid number of vertices
      assert(len(result) >= 3)
      assert(len(result) >= min([len(subject), 4]))

#            for v in result:
#                assert(point_in_polygon(v, clip))

    print 'OK'

#    def tst_non_intersecting(self):
#        from dials.algorithms.polygon import clip
#        from scitbx.array_family import flex

#        for i in range(10000):

#            # Generate nonintersecting polygons
#            subject, target = generate_non_intersecting(target_size=2)
#            rect = ((0, 0), (10, 10))


#            # Do the clipping
#            result = clip.simple_with_rect(
#                flex.vec2_double(subject), rect)

#            print list(result)

#            # Ensure we no vertices
#            assert(len(result) == 0)

#        print 'OK'


class TestTriangleWithTriangle(object):

  def __init__(self):
    pass

  def __call__(self):
    self.tst_intersecting()
    self.tst_non_intersecting()

  def tst_intersecting(self):
    from dials.algorithms.polygon import clip
    for i in range(10000):

      # Generate intersecting polygons
      subject, target = generate_intersecting(3, 3)

      # Do the clipping
      result = clip.triangle_with_triangle(subject, target)

      # Ensure we have roughly valid number of vertices
      assert(len(result) >= 3)
      assert(len(result) >= min([len(subject), len(target)]))

#            for v in result:
#                assert(point_in_polygon(v, clip))

    print 'OK'

  def tst_non_intersecting(self):
    from dials.algorithms.polygon import clip

    for i in range(10000):

      # Generate nonintersecting polygons
      subject, target = generate_non_intersecting(3, 3)

      # Do the clipping
      result = clip.triangle_with_triangle(subject, target)

      # Ensure we no vertices
      assert(len(result) == 0)

    print 'OK'


class TestTriangleWithConvexQuad(object):

  def __init__(self):
    pass

  def __call__(self):
    self.tst_intersecting()
    self.tst_non_intersecting()

  def tst_intersecting(self):
    from dials.algorithms.polygon import clip
    for i in range(10000):

      # Generate intersecting polygons
      subject, target = generate_intersecting(3, 4)

      # Do the clipping
      result = clip.triangle_with_convex_quad(subject, target)

      # Ensure we have roughly valid number of vertices
      assert(len(result) >= 3)
      assert(len(result) >= min([len(subject), len(target)]))

#            for v in result:
#                assert(point_in_polygon(v, clip))

    print 'OK'

  def tst_non_intersecting(self):
    from dials.algorithms.polygon import clip

    for i in range(10000):

      # Generate nonintersecting polygons
      subject, target = generate_non_intersecting(3, 4)

      # Do the clipping
      result = clip.triangle_with_convex_quad(subject, target)

      # Ensure we no vertices
      assert(len(result) == 0)

    print 'OK'


class TestQuadWithTriangle(object):

  def __init__(self):
    pass

  def __call__(self):
    self.tst_intersecting()
    self.tst_non_intersecting()

  def tst_intersecting(self):
    from dials.algorithms.polygon import clip
    for i in range(10000):

      # Generate intersecting polygons
      subject, target = generate_intersecting(4, 3)

      # Do the clipping
      result = clip.quad_with_triangle(subject, target)

      # Ensure we have roughly valid number of vertices
      assert(len(result) >= 3)
      assert(len(result) >= min([len(subject), len(target)]))

#            for v in result:
#                assert(point_in_polygon(v, clip))

    print 'OK'

  def tst_non_intersecting(self):
    from dials.algorithms.polygon import clip

    for i in range(10000):

      # Generate nonintersecting polygons
      subject, target = generate_non_intersecting(4, 3)

      # Do the clipping
      result = clip.quad_with_triangle(subject, target)

      # Ensure we no vertices
      assert(len(result) == 0)

    print 'OK'


class TestQuadWithConvexQuad(object):

  def __init__(self):
    pass

  def __call__(self):
    self.tst_intersecting()
    self.tst_non_intersecting()

  def tst_intersecting(self):
    from dials.algorithms.polygon import clip
    for i in range(10000):

      # Generate intersecting polygons
      subject, target = generate_intersecting(4, 4)

      # Do the clipping
      result = clip.quad_with_convex_quad(subject, target)

      # Ensure we have roughly valid number of vertices
      assert(len(result) >= 3)
      assert(len(result) >= min([len(subject), len(target)]))

#            for v in result:
#                assert(point_in_polygon(v, clip))

    print 'OK'

  def tst_non_intersecting(self):
    from dials.algorithms.polygon import clip

    for i in range(10000):

      # Generate nonintersecting polygons
      subject, target = generate_non_intersecting(4, 4)

      # Do the clipping
      result = clip.quad_with_convex_quad(subject, target)

      # Ensure we no vertices
      assert(len(result) == 0)

    print 'OK'


class TestLineWithRect(object):

  def __init__(self):
    self.box = ((-10, -10), (10, 10))

  def __call__(self):
    from random import uniform
    from dials.algorithms.polygon import clip

    for i in range(1000):

      point1 = uniform(-20, 20), uniform(-20, 20)
      point2 = uniform(-20, 20), uniform(-20, 20)
      line = (point1, point2)
      line, status = clip.line_with_rect(line, self.box)
      if self.intersects(point1, point2):
        assert(status)
      else:
        assert(status == False)


    print 'OK'

  def inbetween(self, x, x0, x1):
    x00 = min([x0, x1])
    x11 = max([x0, x1])
    return x >= x00 and x <= x11

  def intersects(self, point1, point2):
    if self.is_outside(point1) and self.is_outside(point2):
      m = (point2[1] - point1[1]) / (point2[0] - point1[0])
      c = (point1[1] - m * point1[0])
      x = self.box[0][0]
      y = m * x + c
      if (y >= self.box[0][1] and y <= self.box[1][1] and
          self.inbetween(y, point1[1], point2[1])):
        return True
      x = self.box[1][0]
      y = m * x + c
      if (y >= self.box[0][1] and y <= self.box[1][1] and
          self.inbetween(y, point1[1], point2[1])):
        return True
      y = self.box[0][1]
      x = (y - c) / m
      if (x >= self.box[0][0] and x <= self.box[1][0] and
          self.inbetween(x, point1[0], point2[0])):
        return True
      y = self.box[1][1]
      x = (y - c) / m
      if (x >= self.box[0][0] and x <= self.box[0][1] and
          self.inbetween(x, point1[0], point2[0])):
        return True
      return False
    else:
      return True

  def is_outside(self, point):
    return (point[0] < self.box[0][0] or point[1] < self.box[0][1] or
            point[0] > self.box[1][0] or point[1] > self.box[1][1])


class Test(object):

  def __init__(self):
    self.tst_simple_with_convex = TestSimpleWithConvex()
    self.tst_triangle_with_triangle = TestTriangleWithTriangle()
    self.tst_triangle_with_convex_quad = TestTriangleWithConvexQuad()
    self.tst_quad_with_triangle = TestQuadWithTriangle()
    self.tst_quad_with_convex_quad = TestQuadWithConvexQuad()
    self.tst_line_with_rect = TestLineWithRect()
    self.tst_simple_with_rect = TestSimpleWithRect()

  def run(self):
    self.tst_simple_with_convex()
    self.tst_triangle_with_triangle()
    self.tst_triangle_with_convex_quad()
    self.tst_quad_with_triangle()
    self.tst_quad_with_convex_quad()
    self.tst_line_with_rect()
    self.tst_simple_with_rect()


if __name__ == '__main__':
  from dials.test import cd_auto
  with cd_auto(__file__):
    test = Test()
    test.run()

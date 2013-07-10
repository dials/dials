
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

    if subject_size == None:
        subject_size = randint(3, 10)

    if target_size == None:
        target_size = 4

    bbox = (0.0, 10.0, 0.0, 10.0)
    subject = generate_polygon(subject_size, bbox)
    clip = generate_polygon(target_size, bbox)

    return subject, clip

def generate_non_intersecting(subject_size=None, target_size=None):
    from random import randint

    if subject_size == None:
        subject_size = randint(3, 10)

    if target_size == None:
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

class TestTriangleWithTriangle(object):

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
        from scitbx.array_family import flex

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
        from scitbx.array_family import flex
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
        from scitbx.array_family import flex

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
        from scitbx.array_family import flex
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
        from scitbx.array_family import flex

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
        from scitbx.array_family import flex
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
        from scitbx.array_family import flex

        for i in range(10000):

            # Generate nonintersecting polygons
            subject, target = generate_non_intersecting(4, 4)

            # Do the clipping
            result = clip.quad_with_convex_quad(subject, target)

            # Ensure we no vertices
            assert(len(result) == 0)

        print 'OK'


class Test(object):

    def __init__(self):
        self.tst_simple_with_convex = TestSimpleWithConvex()
        self.tst_triangle_with_triangle = TestTriangleWithTriangle()
        self.tst_triangle_with_convex_quad = TestTriangleWithConvexQuad()
        self.tst_quad_with_triangle = TestQuadWithTriangle()
        self.tst_quad_with_convex_quad = TestQuadWithConvexQuad()

    def run(self):
        self.tst_simple_with_convex()
        self.tst_triangle_with_triangle()
        self.tst_triangle_with_convex_quad()
        self.tst_quad_with_triangle()
        self.tst_quad_with_convex_quad()


if __name__ == '__main__':
    test = Test()
    test.run()

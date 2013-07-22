
#def is_inside(p, e1, e2):
#    return (e2[0] - e1[0]) * (p[1] - e1[1]) > (e2[1] - e1[1]) * (p[0] - e1[0])

INSIDE = 0
LEFT = (1 << 0)
RIGHT = (1 << 1)
BOTTOM = (1 << 2)
TOP = (1 << 3)
OUTSIDE = (1 << 4)


def is_inside(p, aabb, side):
    if (side == LEFT):
        return p[0] >= aabb[0][0]
    elif (side == RIGHT):
        return p[0] <= aabb[1][0]
    elif (side == BOTTOM):
        return p[1] >= aabb[0][1]
    elif (side == TOP):
        return p[1] <= aabb[1][1]
    else:
        raise "Bad"

def intersection(p1, p2, aabb, side):
    if side == LEFT:
        x = aabb[0][0]
        y = p1[1] + (p2[1] - p1[1]) * (x - p1[0]) / (p2[0] - p1[0])
        return (x, y)
    elif side == RIGHT:
        x = aabb[1][0]
        y = p1[1] + (p2[1] - p1[1]) * (x - p1[0]) / (p2[0] - p1[0])
        return (x, y)
    elif side == BOTTOM:
        y = aabb[0][1]
        x = p1[0] + (p2[0] - p1[0]) * (y - p1[1]) / (p2[1] - p1[1])
        return (x, y)
    elif side == TOP:
        y = aabb[1][1]
        x = p1[0] + (p2[0] - p1[0]) * (y - p1[1]) / (p2[1] - p1[1])
        return (x, y)
    else:
        raise "Bad"

def clip_line(input, aabb, code):

    # Break if input is empty
    if (len(input) == 0):
        return []

    # Loop through all the edges in the subject polygon, again starting
    # with the last for convenience
    output = []
    p1 = input[-1]
    for p2 in input:

      # If the first point is inside the edge, add the intersection. If
      # the second point is inside the edge, add the point
      if (is_inside(p2, aabb, code)):
        if (not is_inside(p1, aabb, code)):
          output.append(intersection(p1, p2, aabb, code))
        output.append(p2);
      elif (is_inside(p1, aabb, code)):
        output.append(intersection(p1, p2, aabb, code))

      # Advance the subject edge
      p1 = p2

    return output

def sutherland_hodgman(subject, aabb):

    output = subject
    output = clip_line(output, aabb, BOTTOM)
    output = clip_line(output, aabb, RIGHT)
    output = clip_line(output, aabb, TOP)
    output = clip_line(output, aabb, LEFT)

    # Return the clipped polygon vertices
    return output








def display(result, quad, aabb):
    from matplotlib import pylab
    bx = (aabb[0][0], aabb[1][0], aabb[1][0], aabb[0][0])
    by = (aabb[0][1], aabb[0][1], aabb[1][1], aabb[1][1])
    qx, qy = zip(*quad)
    rx, ry = zip(*result)
    bx = bx + (bx[0],)
    by = by + (by[0],)
    qx = qx + (qx[0],)
    qy = qy + (qy[0],)
    rx = rx + (rx[0],)
    ry = ry + (ry[0],)
    pylab.plot(bx, by, color='b')
    pylab.plot(qx, qy, color='y')
    pylab.plot(rx, ry, color='r')
    pylab.show()




aabb = ((0, 0), (10, 10))
##quad = ((3, 3), (5, 3), (5, 5), (3, 5))
quad = ((0, -15), (15, 7.5), (7.5, 15), (-5, 5))
quad = ((5, -2), (12, 5), (5, 12), (-2, 5))


from dials.algorithms.polygon import clip
from scitbx.array_family import flex
from time import time

bx = (aabb[0][0], aabb[1][0], aabb[1][0], aabb[0][0])
by = (aabb[0][1], aabb[0][1], aabb[1][1], aabb[1][1])
convex = zip(bx, by)
poly1 = flex.vec2_double(quad)
poly2 = flex.vec2_double(convex)
rect = aabb

st = time()
for i in range(10000):
    result1 = clip.simple_with_convex(poly1, poly2)
print time() - st

st = time()
for i in range(10000):
#    result2 = sutherland_hodgman(quad, aabb)
    result2 = clip.simple_with_rect(poly1, rect)
print time() - st

print list(result1)
print list(result2)

#poly = sutherland_hodgman(quad, aabb)
#print poly
#display(poly, quad, aabb)

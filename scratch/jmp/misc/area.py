

INSIDE = 0
LEFT = (1 << 0)
RIGHT = (1 << 1)
BOTTOM = (1 << 2)
TOP = (1 << 3)
OUTSIDE = (1 << 4)

def outcode(p, b):
    code = INSIDE
    if (p[0] < b[0][0]): code |= LEFT
    elif(p[0] > b[1][0]): code |= RIGHT
    if (p[1] < b[0][1]): code |= BOTTOM
    elif(p[1] > b[1][1]): code |= TOP
    return code

#def quad_with_aabb_area(poly, aabb):
#
#    p0 = poly[-1]
#    code0 = outcode(p0, aabb)
#    area = 0.0
#
#    for p1 in poly:
#        code1 = outcode(p1, aabb)
#
#        # If both points are in the same zone then add the area
#        # contribution from the whole line between the points
#        if code0 == code1:
#            if code0 != INSIDE:
#                area += p0[0]*p1[1] - p1[0]*p0[1]
#        else:
#            # Get the intersection point between the next zone
#            pass
#
#        # Step forward
#        p0 = p1
#        code0 = code1

#    # Return the intersection area
#    return 0.5 * area

def area_component(i, poly, aabb):

    p0 = poly[i]
    code0 = outcode(p0, aabb)

    area = 0.0
    for j in range(i, len(poly)):
        p1 = poly[j]
        code1 = outcode(p1, aabb)

        if code0 == code1:
            if code0 != INSIDE:
                area += p0[0] * p1[i] - p1[0] * p0[1]


    return area


def quad_with_aabb_area(poly, aabb):

    area = 0.0

    i = 0
    while i < len(poly):

        area, i += area_component(i, poly, aabb)

    # Return the intersection area
    return 0.5 * area


#aabb = ((-10, -10), (10, 10))
#quad = ((0, -15), (15, 7.5), (7.5, 15), (-5, 5))

#from matplotlib import pylab
#bx = [aabb[0][0], aabb[1][0], aabb[1][0], aabb[0][0]]
#by = [aabb[0][1], aabb[0][1], aabb[1][1], aabb[1][1]]
#px, py = zip(*quad)
#pylab.plot(bx, by)
#pylab.plot(px, py)
#pylab.show()

#p0, p1, p2, p3 = quad

#m = (p1[1] - p0[1]) / (p1[0] - p0[0])
#c = p0[1] - m * p0[0]
#y = aabb[0][1]
#x = (y - c) / m
#a0 = 0.5 * x * (p0[1] - y)
#print x, y, a0

#m = (p1[1] - p0[1]) / (p1[0] - p0[0])
#c = p0[1] - m * p0[0]
#x = aabb[1][0]
#y = m * x + c
#a1 = 0.5 * (p1[0] - x) * y
#print x, y, a1

#m = (p2[1] - p1[1]) / (p2[0] - p1[0])
#c = p1[1] - m * p1[0]
#x = aabb[1][0]
#y = m * x + c
#a2 = 0.5 * (p2[0] - x) * y
#print x, y, a2

#m = (p3[1] - p2[1]) / (p3[0] - p2[0])
#c = p2[1] - m * p2[0]
#y = aabb[1][1]
#x = (y - c) / m
#a2 = 0.5 * (y * (p2[0] - aabb[1][0]) - (x - aabb[1][0]) * p2[1])
#print m, c, x, y, a2


aabb = ((0, 0), (10, 10))
#quad = ((5, -10), (12, -8), (-5, -1), (-10, -3))
#quad = ((0, -15), (15, 7.5), (7.5, 15), (-5, 5))
quad = ((3, 3), (5, 3), (5, 5), (3, 5))

from dials.algorithms.polygon import simple_area
from scitbx.array_family import flex
print quad_with_aabb_area(quad, aabb), simple_area(flex.vec2_double(quad))

#p0, p1, p2, p3 = quad

#m0 = (p1[1] - p0[1]) / (p1[0] - p0[0])
#m1 = (p2[1] - p1[1]) / (p2[0] - p1[0])
#m2 = (p3[1] - p2[1]) / (p3[0] - p2[0])
#m3 = (p0[1] - p3[1]) / (p0[0] - p3[0])
#c0 = p0[1] - m0 * p0[0]
#c1 = p1[1] - m1 * p1[0]
#c2 = p2[1] - m2 * p2[0]
#c3 = p3[1] - m3 * p3[0]

#seg = [
# (p0, (aabb[1][0], m0 * aabb[1][0] + c0), (aabb[1][0], aabb[0][1]), (p0[0], aabb[0][1])),
# ((aabb[1][0], m0 * aabb[1][0] + c0), p1, (aabb[1][0], m1 * aabb[1][0] + c1)),
# ((aabb[1][0], m1 * aabb[1][0] + c1), (aabb[0][0], m1 * aabb[0][0] + c1), (aabb[0][0], aabb[0][1]), (aabb[1][0], aabb[0][1])),
# ((aabb[0][0], m1 * aabb[0][0] + c1), p2, p3, (aabb[0][0], m3 * aabb[0][0] + c3)),
# ((aabb[0][0], m3 * aabb[0][0] + c3), p0, (p0[0], aabb[0][1]), (aabb[0][0], aabb[0][1]))
#]

#from dials.algorithms.polygon import simple_area
#from scitbx.array_family import flex
#from matplotlib import pylab
#a = 0
#for s in seg:
#    poly = flex.vec2_double(s)
#    a += simple_area(poly)
#print a

#print simple_area(flex.vec2_double(quad))
##    x, y = zip(*s)
##    x = x + (x[0],)
##    y = y + (y[0],)
##    pylab.plot(x, y)
##pylab.show()


#quad =

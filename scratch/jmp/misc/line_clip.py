
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

def get_point(p, b, outcode, m, c):
    code = INSIDE
    if outcode == INSIDE:
        x = p[0]
        y = p[1]
    elif outcode == LEFT:
        x = b[0][0]
        y = m * x + c
        if (y < b[0][1] or b[1][1] < y):
            code = OUTSIDE
    elif outcode == RIGHT:
        x = b[1][0]
        y = m * x + c
        if (y < b[0][1] or b[1][1] < y):
            code = OUTSIDE
    elif outcode == BOTTOM:
        y = b[0][1]
        x = (y - c) / m
        if (x < b[0][0] or b[1][0] < x):
            code = OUTSIDE
    elif outcode == TOP:
        y = b[1][1]
        x = (y - c) / m
        if (x < b[0][0] or b[1][0] < x):
            code = OUTSIDE
    elif outcode == LEFT | BOTTOM:
        x = b[0][0]
        y = m * x + c
        if y > b[1][1]:
            code = OUTSIDE
        elif y < b[0][1]:
            y = b[0][1]
            x = (y - c) / m
            if (b[1][0] < x):
                code = OUTSIDE
    elif outcode == LEFT | TOP:
        x = b[0][0]
        y = m * x + c
        if y < b[0][1]:
            code = OUTSIDE
        elif y > b[1][1]:
            y = b[1][1]
            x = (y - c) / m
            if (b[1][0] < x):
                code = OUTSIDE
    elif outcode == RIGHT | BOTTOM:
        x = b[1][0]
        y = m * x + c
        if y > b[1][1]:
            code = OUTSIDE
        elif y < b[0][1]:
            y = b[0][1]
            x = (y - c) / m
            if (b[0][0] > x):
                code = OUTSIDE
    elif outcode == RIGHT | TOP:
        x = b[1][0]
        y = m * x + c
        if y < b[0][1]:
            code = OUTSIDE
        elif y > b[1][1]:
            y = b[1][1]
            x = (y - c) / m
            if (b[0][0] > x):
                code = OUTSIDE
    else:
        raise RuntimeError("Noooo")

    return (x, y), code

def clip_line(box, line):

    outcode1 = outcode(box[0], box[1], box[2], box[3], line[0], line[1])
    outcode2 = outcode(box[0], box[1], box[2], box[3], line[2], line[3])

    if not (outcode1 or outcode2):
        return line
    if outcode1 & outcode2:
        return None

    m = (line[3] - line[1]) / (line[2] - line[0])
    c = (line[1] - line[0] * m)

    point1, outcode1 = get_point(box, (line[0], line[1]), outcode1, m, c)
    point2, outcode2 = get_point(box, (line[2], line[3]), outcode2, m, c)

    if not (outcode1 or outcode2):
        return point1 + point2
    else:
        return None

def get_corner(code):
    corner = 0
    if code == LEFT | BOTTOM:
        corner = 1
    elif code == LEFT | TOP:
        corner = 2
    elif code == RIGHT | BOTTOM:
        corner = 3
    elif code == RIGHT | TOP:
        corner = 4
    return corner

LB = (1 << 0)
RB = (1 << 1)
LT = (1 << 2)
RT = (1 << 3)


def add_corner(code0, code1):

    corner = 0
    if code0 & LEFT:
        if not (code1 & LEFT):
            if code1 & TOP:
                corner |= LT
            elif code1 & BOTTOM:
                corner |= LB
    if code0 & RIGHT:
        if not (code1 & RIGHT):
            if code1 & TOP:
                corner |= RT
            elif code1 & BOTTOM:
                corner |= RB
    if code0 & TOP:
        if not (code1 & TOP):
            if code1 & LEFT:
                corner |= LT
            elif code1 & RIGHT:
                corner |= RT
    if code0 & BOTTOM:
        if not (code1 & BOTTOM):
            if code1 & LEFT:
                corner |= LB
            elif code1 & RIGHT:
                corner |= RB

    return corner

def convex_quad_with_aabb(quad, aabb):

    p0, p1, p2, p3 = quad
    b0 = aabb[0][0], aabb[0][1]
    b1 = aabb[1][0], aabb[0][1]
    b2 = aabb[1][0], aabb[1][1]
    b3 = aabb[0][0], aabb[1][1]

    code0 = outcode(p0, aabb)
    code1 = outcode(p1, aabb)
    code2 = outcode(p2, aabb)
    code3 = outcode(p3, aabb)

    m0 = (p1[1] - p0[1]) / (p1[0] - p0[0])
    m1 = (p2[1] - p1[1]) / (p2[0] - p1[0])
    m2 = (p3[1] - p2[1]) / (p3[0] - p2[0])
    m3 = (p0[1] - p3[1]) / (p0[0] - p3[0])
    c0 = p0[1] - m0 * p0[0]
    c1 = p1[1] - m1 * p1[0]
    c2 = p2[1] - m2 * p2[0]
    c3 = p3[1] - m3 * p3[0]

    p00, cc00 = get_point(p0, aabb, code0, m0, c0)
    p01, cc01 = get_point(p1, aabb, code1, m0, c0)
    p11, cc11 = get_point(p1, aabb, code1, m1, c1)
    p12, cc12 = get_point(p2, aabb, code2, m1, c1)
    p22, cc22 = get_point(p2, aabb, code2, m2, c2)
    p23, cc23 = get_point(p3, aabb, code3, m2, c2)
    p33, cc33 = get_point(p3, aabb, code3, m3, c3)
    p30, cc30 = get_point(p0, aabb, code0, m3, c3)

    poly = []

    bx = (aabb[0][0], aabb[1][0], aabb[1][0], aabb[0][0])
    by = (aabb[0][1], aabb[0][1], aabb[1][1], aabb[1][1])
    pc = zip(bx, by)

    if not (code0 & code1) and not (cc00 or cc01):
        poly.append(p00)
        if code1 != INSIDE:
            poly.append(p01)
    else:
        corner = add_corner(code0, code1)
        if corner & LB:
            poly.append(pc[0])
        if corner & RB:
            poly.append(pc[1])
        if corner & RT:
            poly.append(pc[2])
        if corner & LT:
            poly.append(pc[3])

    if not (code1 & code2) and not (cc11 or cc12):
        poly.append(p11)
        if code2 != INSIDE:
            poly.append(p12)
    else:
        corner = add_corner(code1, code2)
        if corner & LB:
            poly.append(pc[0])
        if corner & RB:
            poly.append(pc[1])
        if corner & RT:
            poly.append(pc[2])
        if corner & LT:
            poly.append(pc[3])

    if not (code2 & code3) and not (cc22 or cc23):
        poly.append(p22)
        if code3 != INSIDE:
            poly.append(p23)
    else:
        corner = add_corner(code2, code3)
        if corner & LB:
            poly.append(pc[0])
        if corner & RB:
            poly.append(pc[1])
        if corner & RT:
            poly.append(pc[2])
        if corner & LT:
            poly.append(pc[3])

    if not (code3 & code0) and not (cc33 or cc30):
        poly.append(p33)
        if code0 != INSIDE:
            poly.append(p30)
    else:
        corner = add_corner(code3, code0)
        if corner & LB:
            poly.append(pc[0])
        if corner & RB:
            poly.append(pc[1])
        if corner & RT:
            poly.append(pc[2])
        if corner & LT:
            poly.append(pc[3])

    return poly




def add_corner_point(code0, code1):

    corner = 0
    if code0 & LEFT:
        if not (code1 & LEFT):
            if code1 & TOP:
                corner |= LT
            elif code1 & BOTTOM:
                corner |= LB
    if code0 & RIGHT:
        if not (code1 & RIGHT):
            if code1 & TOP:
                corner |= RT
            elif code1 & BOTTOM:
                corner |= RB
    if code0 & TOP:
        if not (code1 & TOP):
            if code1 & LEFT:
                corner |= LT
            elif code1 & RIGHT:
                corner |= RT
    if code0 & BOTTOM:
        if not (code1 & BOTTOM):
            if code1 & LEFT:
                corner |= LB
            elif code1 & RIGHT:
                corner |= RB

    return corner

def convex_quad_with_aabb_area(quad, aabb):

    p0, p1, p2, p3 = quad
    b0 = aabb[0][0], aabb[0][1]
    b1 = aabb[1][0], aabb[0][1]
    b2 = aabb[1][0], aabb[1][1]
    b3 = aabb[0][0], aabb[1][1]

    code0 = outcode(p0, aabb)
    code1 = outcode(p1, aabb)
    code2 = outcode(p2, aabb)
    code3 = outcode(p3, aabb)

    m0 = (p1[1] - p0[1]) / (p1[0] - p0[0])
    m1 = (p2[1] - p1[1]) / (p2[0] - p1[0])
    m2 = (p3[1] - p2[1]) / (p3[0] - p2[0])
    m3 = (p0[1] - p3[1]) / (p0[0] - p3[0])
    c0 = p0[1] - m0 * p0[0]
    c1 = p1[1] - m1 * p1[0]
    c2 = p2[1] - m2 * p2[0]
    c3 = p3[1] - m3 * p3[0]

    p00, cc00 = get_point(p0, aabb, code0, m0, c0)
    p01, cc01 = get_point(p1, aabb, code1, m0, c0)
    p11, cc11 = get_point(p1, aabb, code1, m1, c1)
    p12, cc12 = get_point(p2, aabb, code2, m1, c1)
    p22, cc22 = get_point(p2, aabb, code2, m2, c2)
    p23, cc23 = get_point(p3, aabb, code3, m2, c2)
    p33, cc33 = get_point(p3, aabb, code3, m3, c3)
    p30, cc30 = get_point(p0, aabb, code0, m3, c3)

    area = 0.0

    bx = (aabb[0][0], aabb[1][0], aabb[1][0], aabb[0][0])
    by = (aabb[0][1], aabb[0][1], aabb[1][1], aabb[1][1])
    pc = zip(bx, by)

    if not (code0 & code1) and not (cc00 or cc01):
        area += p00[0] * p01[1] - p01[0] * p00[1]
        pl = p01
    else:
        corner = add_corner_point(code0, code1)
        if corner & LB:
            poly.append(pc[0])
        if corner & RB:
            poly.append(pc[1])
        if corner & RT:
            poly.append(pc[2])
        if corner & LT:
            poly.append(pc[3])

    if not (code1 & code2) and not (cc11 or cc12):
        area += p11[0] * p12[1] - p12[0] * p11[1]
        pl = p12
    else:
        corner = add_corner_point(code1, code2)
        if corner & LB:
            poly.append(pc[0])
        if corner & RB:
            poly.append(pc[1])
        if corner & RT:
            poly.append(pc[2])
        if corner & LT:
            poly.append(pc[3])

    if not (code2 & code3) and not (cc22 or cc23):
        area += p22[0] * p23[1] - p23[0] * p22[1]
        pl = p23
    else:
        corner = add_corner_point(code2, code3)
        if corner & LB:
            poly.append(pc[0])
        if corner & RB:
            poly.append(pc[1])
        if corner & RT:
            poly.append(pc[2])
        if corner & LT:
            poly.append(pc[3])

    if not (code3 & code0) and not (cc33 or cc30):
        area += p33[0] * p30[1] - p30[0] * p33[1]
        pl = p30
    else:
        corner = add_corner_point(code3, code0)
        if corner & LB:
            poly.append(pc[0])
        if corner & RB:
            poly.append(pc[1])
        if corner & RT:
            poly.append(pc[2])
        if corner & LT:
            poly.append(pc[3])

    return 0.5 * area






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

class Test(object):

    def run(self):

        self.tst_single_zone()
        self.tst_outside_spanning()

    def tst_single_zone(self):

        # Box
        aabb = ((0, 0), (10, 10))

        # Inside
        quad = ((5, 2.5), (7.5, 5), (5, 7.5), (2.5, 5))
        poly = convex_quad_with_aabb(quad, aabb)
        assert(len(poly) == 4)

        # Left
        quad = ((-10, 2.5), (-7.5, 5), (-10, 7.5), (-12.5, 5))
        poly = convex_quad_with_aabb(quad, aabb)
        assert(len(poly) == 0)

        # Right
        quad = ((20, 2.5), (17.5, 5), (20, 7.5), (22.5, 5))
        poly = convex_quad_with_aabb(quad, aabb)
        assert(len(poly) == 0)

        # TOP
        quad = ((5, 12.5), (7.5, 15), (5, 17.5), (2.5, 15))
        poly = convex_quad_with_aabb(quad, aabb)
        assert(len(poly) == 0)

        # BOTTOM
        quad = ((5, -2.5), (7.5, -5), (5, -7.5), (2.5, -5))
        poly = convex_quad_with_aabb(quad, aabb)
        assert(len(poly) == 0)

        # Bottom left
        quad = ((-5, -2.5), (-7.5, -5), (-5, -7.5), (-2.5, -5))
        poly = convex_quad_with_aabb(quad, aabb)
        assert(len(poly) == 0)

        # Bottom Right
        quad = ((15, -2.5), (17.5, -5), (15, -7.5), (12.5, -5))
        poly = convex_quad_with_aabb(quad, aabb)
        assert(len(poly) == 0)

        # Top left
        quad = ((-5, 12.5), (-7.5, 15), (-5, 17.5), (-2.5, 15))
        poly = convex_quad_with_aabb(quad, aabb)
        assert(len(poly) == 0)

        # Top Right
        quad = ((15, 12.5), (17.5, 15), (15, 17.5), (12.5, 15))
        poly = convex_quad_with_aabb(quad, aabb)
        assert(len(poly) == 0)

        # Test passed
        print 'OK'

    def tst_outside_spanning(self):

        # Box
        aabb = ((0, 0), (10, 10))

        # Bottom left to bottom right
        quad = ((5, -2.5), (-7.5, -5), (5, -7.5), (17.5, -5))
        poly = convex_quad_with_aabb(quad, aabb)
        assert(len(poly) == 2)

        # Top left to Top right
        quad = ((5, 12.5), (-7.5, 15), (5, 17.5), (17.5, 15))
        poly = convex_quad_with_aabb(quad, aabb)
        assert(len(poly) == 2)

        # Test passed
        print 'OK'

test = Test()
test.run()


#
#
#aabb = ((0, 0), (10, 10))
##quad = ((3, 3), (5, 3), (5, 5), (3, 5))
#quad = ((0, -15), (15, 7.5), (7.5, 15), (-5, 5))

#result = convex_quad_with_aabb(quad, aabb)
#print result
#from matplotlib import pylab
#bx = (aabb[0][0], aabb[1][0], aabb[1][0], aabb[0][0])
#by = (aabb[0][1], aabb[0][1], aabb[1][1], aabb[1][1])
#qx, qy = zip(*quad)
#rx, ry = zip(*result)
#bx = bx + (bx[0],)
#by = by + (by[0],)
#qx = qx + (qx[0],)
#qy = qy + (qy[0],)
#rx = rx + (rx[0],)
#ry = ry + (ry[0],)
#pylab.plot(bx, by, color='b')
#pylab.plot(qx, qy, color='y')
#pylab.plot(rx, ry, color='r')
#pylab.show()

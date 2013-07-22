
def is_inside(p, e1, e2):
    return (e2[0] - e1[0]) * (p[1] - e1[1]) > (e2[1] - e1[1]) * (p[0] - e1[0])

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


def intersection(p1, p2, e1, e2):
    dc = (e1[0] - e2[0], e1[1] - e2[1])
    dp = (p1[0] - p2[0], p1[1] - p2[1])
    n1 = e1[0] * e2[1] - e1[1] * e2[0]
    n2 = p1[0] * p2[1] - p1[1] * p2[0]
    n3 = (dc[0] * dp[1] - dc[1] * dp[0])
    assert(n3 != 0.0)
    n4 = 1.0 / n3
    return ((n1 * dp[0] - n2 * dc[0]) * n4, (n1 * dp[1] - n2 * dc[1]) * n4)


def clip_line(input, e1, e2, aabb, code):

    # Break if input is empty
    if (len(input) == 0):
        return []

    # Loop through all the edges in the subject polygon, again starting
    # with the last for convenience
    output = []
    p1 = input[-1]
    code1 = outcode(p1, aabb)
    for p2 in input:
      code2 = outcode(p2, aabb)


      # If the first point is inside the edge, add the intersection. If
      # the second point is inside the edge, add the point
      if (is_inside(p2, e1, e2)):
        if (not is_inside(p1, e1, e2)):
          output.append(intersection(p1, p2, e1, e2))
        output.append(p2);
      elif (is_inside(p1, e1, e2)):
        output.append(intersection(p1, p2, e1, e2))

      # Advance the subject edge
      p1 = p2

    return output

def sutherland_hodgman(subject, aabb):

    tx = aabb[0][0], aabb[1][0], aabb[1][0], aabb[0][0]
    ty = aabb[0][1], aabb[0][1], aabb[1][1], aabb[1][1]
    target = zip(tx, ty)

    # Ensure polygons are valid
    assert(len(subject) >= 3 and len(target) >= 3)

    output = subject


    output = clip_line(output, target[0], target[1], aabb, 0)
    output = clip_line(output, target[1], target[2], aabb, 0)
    output = clip_line(output, target[2], target[3], aabb, 0)
    output = clip_line(output, target[3], target[0], aabb, 0)


#    # Loop through all the edges of the clip polygon, starting with the
#    # last edge for conveneince
#    e1 = target[-1]
#    for e2 in target:

#      # Break if input is empty
#      if (len(output) == 0):
#          break

#      output = clip_line(output, e1, e2, aabb, 0)
#
#      # Advance the target edge
#      e1 = e2

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

poly = sutherland_hodgman(quad, aabb)
display(poly, quad, aabb)

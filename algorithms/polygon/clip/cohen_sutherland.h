/*
 * cohen_sutherland.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_POLYGON_CLIP_COHEN_SUTHERLAND_H
#define DIALS_ALGORITHMS_POLYGON_CLIP_COHEN_SUTHERLAND_H

#include <scitbx/array_family/tiny_types.h>
#include <scitbx/array_family/small.h>
#include <dials/error.h>

namespace dials { namespace algorithms { namespace polygon { namespace clip {

  using scitbx::af::double2;
  using scitbx::af::double4;
  using scitbx::af::small;

  // Constants used in algorithm
  enum {
    INSIDE = 0,
    LEFT = (1 << 0),
    RIGHT = (1 << 1),
    BOTTOM = (1 << 2),
    TOP = (1 << 3),
    OUTSIDE = (1 << 4),
  };

  /**
   * Get the outcode for the point in the box
   * @param p The point
   * @param b The box
   * @returns The code describing the location of the point relative to the box.
   */
  template <typename PointType, typename BoxType>
  int cohen_sutherland_outcode(const PointType &p, const BoxType &b) {
    int code = INSIDE;
    if (p[0] < b[0][0])
      code |= LEFT;
    else if (p[0] > b[1][0])
      code |= RIGHT;
    if (p[1] < b[0][1])
      code |= BOTTOM;
    else if (p[1] > b[1][1])
      code |= TOP;
    return code;
  }

  /**
   * A large switch statement to quickly find the intersection point if
   * any of a point with the box.
   * @param p The point
   * @param aabb The box
   * @param code The input code
   * @param m The line gradient
   * @param c The line intersect with 0
   * @retuns The return code.
   */
  template <typename PointType, typename BoxType>
  int cohen_sutherland_intersection(PointType &p,
                                    const BoxType &aabb,
                                    int code,
                                    double m,
                                    double c) {
    int retcode = INSIDE;
    double x = 0, y = 0;

    // Handle the point and get the intersection
    switch (code) {
    case INSIDE:
      break;
    case LEFT:
      x = aabb[0][0];
      y = m * x + c;
      if (y < aabb[0][1] || aabb[1][1] < y) retcode = OUTSIDE;
      break;
    case RIGHT:
      x = aabb[1][0];
      y = m * x + c;
      if (y < aabb[0][1] || aabb[1][1] < y) retcode = OUTSIDE;
      break;
    case BOTTOM:
      y = aabb[0][1];
      x = (y - c) / m;
      if (x < aabb[0][0] || aabb[1][0] < x) retcode = OUTSIDE;
      break;
    case TOP:
      y = aabb[1][1];
      x = (y - c) / m;
      if (x < aabb[0][0] || aabb[1][0] < x) retcode = OUTSIDE;
      break;
    case LEFT | BOTTOM:
      x = aabb[0][0];
      y = m * x + c;
      if (y > aabb[1][1]) {
        retcode = OUTSIDE;
      } else if (y < aabb[0][1]) {
        y = aabb[0][1];
        x = (y - c) / m;
        if (aabb[1][0] < x) retcode = OUTSIDE;
      }
      break;
    case LEFT | TOP:
      x = aabb[0][0];
      y = m * x + c;
      if (y < aabb[0][1])
        retcode = OUTSIDE;
      else if (y > aabb[1][1]) {
        y = aabb[1][1];
        x = (y - c) / m;
        if (aabb[1][0] < x) retcode = OUTSIDE;
      }
      break;
    case RIGHT | BOTTOM:
      x = aabb[1][0];
      y = m * x + c;
      if (y > aabb[1][1]) {
        retcode = OUTSIDE;
      } else if (y < aabb[0][1]) {
        y = aabb[0][1];
        x = (y - c) / m;
        if (aabb[0][0] > x) retcode = OUTSIDE;
      }
      break;
    case RIGHT | TOP:
      x = aabb[1][0];
      y = m * x + c;
      if (y < aabb[0][1]) {
        retcode = OUTSIDE;
      } else if (y > aabb[1][1]) {
        y = aabb[1][1];
        x = (y - c) / m;
        if (aabb[0][0] > x) retcode = OUTSIDE;
      }
      break;
    default:
      throw DIALS_ERROR("Shouldn't reach this point!");
      break;
    };

    // Return the x, y coordinates
    p[0] = x;
    p[1] = y;

    // Return the code
    return retcode;
  }

  /**
   * Do the cohen sutherland algorithm for finding the intersection of a line
   * with an axis aligned bounding box.
   * @param line The line to intersect
   * @param aabb The axis aligned bounding box
   * @returns The intersecting line
   */
  template <typename LineType, typename BoxType>
  std::pair<LineType, bool> cohen_sutherland_single(const LineType &line,
                                                    const BoxType &aabb) {
    typedef typename LineType::value_type point_type;

    // Get the points on the line
    point_type p1(line[0][0], line[0][1]);
    point_type p2(line[1][0], line[1][1]);

    // Get the outcodes for the two points on the line segment
    int code1 = cohen_sutherland_outcode(p1, aabb);
    int code2 = cohen_sutherland_outcode(p2, aabb);

    // If both ends are inside the box then return the line
    if (!(code1 || code2)) {
      return std::make_pair(LineType(p1, p2), true);
    }

    // If the line can't possible intersect then return Null
    if (code1 & code2) {
      return std::make_pair(LineType(), false);
    }

    // Get the line coefficients
    double m = (p2[1] - p1[1]) / (p2[0] - p1[0]);
    double c = (p1[1] - p1[0] * m);

    // Get the intersection of the line
    code1 = cohen_sutherland_intersection(p1, aabb, code1, m, c);
    code2 = cohen_sutherland_intersection(p2, aabb, code2, m, c);

    // Check if the line intersects and return
    if (!(code1 || code2)) {
      return std::make_pair(LineType(p1, p2), true);
    }

    // Otherwise return Null
    return std::make_pair(LineType(), false);
  }

}}}}  // namespace dials::algorithms::polygon::clip

#endif /* DIALS_ALGORITHMS_POLYGON_CLIP_COHEN_SUTHERLAND_H */

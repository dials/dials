/*
 * sutherland_hodgman.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_POLYGON_CLIPPING_SUTHERLAND_HODGMAN_H
#define DIALS_ALGORITHMS_POLYGON_CLIPPING_SUTHERLAND_HODGMAN_H

#include <dials/error.h>

namespace dials { namespace algorithms { namespace polygon {

  /**
   * Check if a point is "inside" and edge. A point is defined as inside
   * if it lies on the same side of the edge as the remainder of the polygon.
   * Since we assume that the polygon is given in anticlockwise order, this is
   * equivalent to testing whether the point lies to the left of the line.
   * @param p The point
   * @param e1 The first edge point
   * @param e2 The second edge point
   */
  template <typename T>
  bool is_inside(const T &p, const T &e1, const T &e2) {
    return (e2[0] - e1[0]) * (p[1] - e1[1]) > (e2[1] - e1[1]) * (p[0] - e1[0]);
  }

  /**
   * Return the intersection of a line segment with an infinite edge.
   * @param p1 The first point
   * @param p2 The second point
   * @param e1 The first edge point
   * @param e2 The second edge point
   * @returns The intersection
   */
  template <typename T>
  T intersection(const T &p1, const T &p2, const T &e1, const T &e2) {
    T dc(e1[0] - e2[0], e1[1] - e2[1]);
    T dp(p1[0] - p2[0], p1[1] - p2[1]);
    double n1 = e1[0] * e2[1] - e1[1] * e2[0];
    double n2 = p1[0] * p2[1] - p1[1] * p2[0];
    double n3 = (dc[0] * dp[1] - dc[1] * dp[0]);
    DIALS_ASSERT(n3 != 0.0);
    double n4 = 1.0 / n3;
    return T((n1 * dp[0] - n2 * dc[0]) * n4, (n1 * dp[1] - n2 * dc[1]) * n4);
  }

  /**
   * Clip the subject polygon by the clip polygon and return the resulting
   * polygon which will have at most the sum of the edges of the inputs.
   * @param subject The subject polygon
   * @param target The clip polygon
   * @returns The resulting clipped polygon.
   */
  template <typename PolygonType>
  PolygonType sutherland_hodgman(const PolygonType &subject,
                                 const PolygonType &target) {
    // Get the polygon and vertex type
    typedef PolygonType polygon_type;
    typedef typename PolygonType::value_type vertex_type;

    // Ensure polygons are valid
    DIALS_ASSERT(subject.size() >= 3 && target.size() >= 3);

    // Create the buffers
    polygon_type input;
    polygon_type output(subject.size());

    // Copy the subject vertices
    for (std::size_t i = 0; i < output.size(); ++i) {
      output[i] = subject[i];
    }

    // Loop through all the edges of the clip polygon, starting with the
    // last edge for conveneince
    vertex_type e1 = target[target.size() - 1];
    for (std::size_t j = 0; j < target.size(); ++j) {
      vertex_type e2 = target[j];

      // Swap the buffers and vertex counts
      std::swap(input, output);
      output.clear();

      // Break if input is empty
      if (input.size() == 0) {
        break;
      }

      // Loop through all the edges in the subject polygon, again starting
      // with the last for convenience
      vertex_type p1 = input[input.size() - 1];
      for (std::size_t i = 0; i < input.size(); ++i) {
        vertex_type p2 = input[i];

        // If the first point is inside the edge, add the intersection. If
        // the second point is inside the edge, add the point
        if (is_inside(p2, e1, e2)) {
          if (!is_inside(p1, e1, e2)) {
            output.push_back(intersection(p1, p2, e1, e2));
          }
          output.push_back(p2);
        } else if (is_inside(p1, e1, e2)) {
          output.push_back(intersection(p1, p2, e1, e2));
        }

        // Advance the subject edge
        p1 = p2;
      }

      // Advance the target edge
      e1 = e2;
    }

    // Return the clipped polygon vertices
    return output;
  }

  /**
   * Clip the subject polygon by the clip polygon and return the resulting
   * polygon which will have at most the sum of the edges of the inputs.
   * @param subject The subject polygon
   * @param target The clip polygon
   * @returns The resulting clipped polygon.
   */
  template <typename SubjectPolygonType,
            typename TargetPolygonType,
            typename OutputPolygonType,
            int MAX_VERTICES>
  OutputPolygonType sutherland_hodgman_simple_convex(const SubjectPolygonType &subject,
                                                     const TargetPolygonType &target) {
    // Get the polygon and vertex type
    typedef typename OutputPolygonType::value_type vertex_type;

    // Create the buffers
    OutputPolygonType input(MAX_VERTICES);
    OutputPolygonType output(MAX_VERTICES);

    // Copy the subject vertices to the output buffer
    std::size_t input_size = 0, output_size = subject.size();
    for (std::size_t i = 0; i < subject.size(); ++i) {
      output[i] = subject[i];
    }

    // Loop through all the edges of the clip polygon, starting with the
    // last edge for conveneince
    vertex_type e1 = target[target.size() - 1];
    for (std::size_t j = 0; j < target.size(); ++j) {
      vertex_type e2 = target[j];

      // Swap the buffers and vertex counts
      std::swap(input, output);
      input_size = output_size;
      output_size = 0;

      // Break if input is empty
      if (input_size == 0) {
        break;
      }

      // Loop through all the edges in the subject polygon, again starting
      // with the last for convenience
      vertex_type p1 = input[input_size - 1];
      for (std::size_t i = 0; i < input_size; ++i) {
        vertex_type p2 = input[i];

        // If the first point is inside the edge, add the intersection. If
        // the second point is inside the edge, add the point
        if (is_inside(p2, e1, e2)) {
          if (!is_inside(p1, e1, e2)) {
            output[output_size++] = intersection(p1, p2, e1, e2);
          }
          output[output_size++] = p2;
        } else if (is_inside(p1, e1, e2)) {
          output[output_size++] = intersection(p1, p2, e1, e2);
        }

        // Advance the subject edge
        p1 = p2;
      }

      // Advance the clip edge
      e1 = e2;
    }

    // Return the clipped polygon vertices
    output.resize(output_size);
    return output;
  }

  // Enum of sides
  enum {
    LEFT = (1 << 0),
    RIGHT = (1 << 1),
    BOTTOM = (1 << 2),
    TOP = (1 << 3),
  };

  /**
   * Check that the point is inside or outside the rectangle
   * @param p The point
   * @param r The rectangle
   * @returns True/False inside or not
   */
  template <int side, typename PointType, typename RectType>
  bool is_inside_rect(const PointType &p, const RectType &r) {
    bool inside = false;
    if (side == LEFT) {
      inside = p[0] >= r[0][0];
    } else if (side == RIGHT) {
      inside = p[0] <= r[1][0];
    } else if (side == BOTTOM) {
      inside = p[1] >= r[0][1];
    } else if (side == TOP) {
      inside = p[1] <= r[1][1];
    } else {
      throw DIALS_ERROR("Unreachable");
    }
    return inside;
  }

  /**
   * Get the intersection of a line with the side of a rectangle.
   * @param p1 The first point on the line.
   * @param p2 The second point on the line.
   * @param r The rectangle
   * @returns The intersection point.
   */
  template <int side, typename PointType, typename RectType>
  PointType intersection_rect(const PointType &p1,
                              const PointType &p2,
                              const RectType &r) {
    PointType p;
    if (side == LEFT) {
      p[0] = r[0][0];
      p[1] = p1[1] + (p2[1] - p1[1]) * (p[0] - p1[0]) / (p2[0] - p1[0]);
    } else if (side == RIGHT) {
      p[0] = r[1][0];
      p[1] = p1[1] + (p2[1] - p1[1]) * (p[0] - p1[0]) / (p2[0] - p1[0]);
    } else if (side == BOTTOM) {
      p[1] = r[0][1];
      p[0] = p1[0] + (p2[0] - p1[0]) * (p[1] - p1[1]) / (p2[1] - p1[1]);
    } else if (side == TOP) {
      p[1] = r[1][1];
      p[0] = p1[0] + (p2[0] - p1[0]) * (p[1] - p1[1]) / (p2[1] - p1[1]);
    } else {
      throw DIALS_ERROR("Unreachable");
    }
    return p;
  }

  /**
   * Clip the polygon against a side of the rectangle.
   * @param result The resulting polygon.
   * @param poly The polygon to clip.
   * @param rect The rectangle to clip against
   */
  template <int side, typename PolygonType, typename RectType>
  void sutherland_hodgman_rect_line(PolygonType &result,
                                    const PolygonType &poly,
                                    const RectType &rect) {
    typedef typename PolygonType::value_type PointType;

    if (poly.size() == 0) {
      return;
    }

    // Loop through all the edges in the subject polygon, again starting
    // with the last for convenience
    PointType p1 = poly[poly.size() - 1];
    for (std::size_t i = 0; i < poly.size(); ++i) {
      PointType p2 = poly[i];

      // If the first point is inside the edge, add the intersection. If
      // the second point is inside the edge, add the point
      if (is_inside_rect<side>(p2, rect)) {
        if (!is_inside_rect<side>(p1, rect)) {
          result.push_back(intersection_rect<side>(p1, p2, rect));
        }
        result.push_back(p2);
      } else if (is_inside_rect<side>(p1, rect)) {
        result.push_back(intersection_rect<side>(p1, p2, rect));
      }

      // Advance the subject edge
      p1 = p2;
    }
  }

  /**
   * Clip a simple polygon against a rectangle using the sutherland hodgman
   * algorithm.
   * @param poly The polygon to clip
   * @param rect The rectangle to clip against
   * @returns The clipped polygon.
   */
  template <typename PolygonType, typename RectType>
  PolygonType sutherland_hodgman_rect(const PolygonType &poly, const RectType &rect) {
    // Clip along each rect line
    PolygonType result1(poly.size()), result2(poly.size());
    result1.clear();
    result2.clear();
    sutherland_hodgman_rect_line<BOTTOM>(result1, poly, rect);
    sutherland_hodgman_rect_line<RIGHT>(result2, result1, rect);
    result1.clear();
    sutherland_hodgman_rect_line<TOP>(result1, result2, rect);
    result2.clear();
    sutherland_hodgman_rect_line<LEFT>(result2, result1, rect);

    // Return the clipped polygon vertices
    return result2;
  }

}}}  // namespace dials::algorithms::polygon

#endif /* DIALS_ALGORITHMS_POLYGON_CLIPPING_SUTHERLAND_HODGMAN_H */

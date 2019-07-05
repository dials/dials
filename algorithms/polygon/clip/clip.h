/*
 * clipping.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_POLYGON_CLIPPING_CLIPPING_H
#define DIALS_ALGORITHMS_POLYGON_CLIPPING_CLIPPING_H

#include <scitbx/vec2.h>
#include <scitbx/array_family/tiny.h>
#include <scitbx/array_family/small.h>
#include <dials/algorithms/polygon/clip/cohen_sutherland.h>
#include <dials/algorithms/polygon/clip/sutherland_hodgman.h>
#include <dials/array_family/scitbx_shared_and_versa.h>

namespace dials { namespace algorithms { namespace polygon { namespace clip {

  using scitbx::vec2;
  using scitbx::af::small;
  using scitbx::af::tiny;

  // Convenient typedefs
  typedef tiny<vec2<double>, 2> vert2;
  typedef tiny<vec2<double>, 3> vert3;
  typedef tiny<vec2<double>, 4> vert4;
  typedef small<vec2<double>, 6> vert6;
  typedef small<vec2<double>, 7> vert7;
  typedef small<vec2<double>, 8> vert8;
  typedef af::shared<vec2<double> > shared_vec2_double;

  /**
   * Clip a simple polygon with a convex polygon using the sutherland_hodman
   * algorithm.
   * @param subject The subject polygon
   * @param target The clip polygon
   * @returns The intersecting polygon
   */
  inline shared_vec2_double simple_with_convex(const shared_vec2_double &subject,
                                               const shared_vec2_double &target) {
    return sutherland_hodgman<shared_vec2_double>(subject, target);
  }

  /**
   * Clip a simple polygon with a rectangle using the sutherland_hodman
   * algorithm.
   * @param poly The subject polygon
   * @param rect The clip polygon
   * @returns The intersecting polygon
   */
  inline shared_vec2_double simple_with_rect(const shared_vec2_double &poly,
                                             const tiny<vec2<double>, 2> &rect) {
    return sutherland_hodgman_rect(poly, rect);
  }

  /**
   * Clip a triangle with a triangle using the sutherland_hodman
   * algorithm.
   * @param subject The subject polygon
   * @param target The clip polygon
   * @returns The intersecting polygon
   */
  inline vert6 triangle_with_triangle(const vert3 &subject, const vert3 &target) {
    return sutherland_hodgman_simple_convex<vert3, vert3, vert6, 6>(subject, target);
  }

  /**
   * Clip a triangle with a convex quad using the sutherland_hodman
   * algorithm.
   * @param subject The subject polygon
   * @param target The clip polygon
   * @returns The intersecting polygon
   */
  inline vert7 triangle_with_convex_quad(const vert3 &subject, const vert4 &target) {
    return sutherland_hodgman_simple_convex<vert3, vert4, vert7, 7>(subject, target);
  }

  /**
   * Clip a quad with a triangle using the sutherland_hodman
   * algorithm.
   * @param subject The subject polygon
   * @param target The clip polygon
   * @returns The intersecting polygon
   */
  inline vert7 quad_with_triangle(const vert4 &subject, const vert3 &target) {
    return sutherland_hodgman_simple_convex<vert4, vert3, vert7, 7>(subject, target);
  }

  /**
   * Clip a quad with a convex quad using the sutherland_hodman
   * algorithm.
   * @param subject The subject polygon
   * @param target The clip polygon
   * @returns The intersecting polygon
   */
  inline vert8 quad_with_convex_quad(const vert4 &subject, const vert4 &target) {
    return sutherland_hodgman_simple_convex<vert4, vert4, vert8, 8>(subject, target);
  }

  /**
   * Clip a line with an axis aligned bounding box.
   * @param line The line to clip
   * @param rect The box
   * @returns The intersecting line
   */
  inline std::pair<vert2, bool> line_with_rect(const vert2 &line, const vert2 &rect) {
    return cohen_sutherland_single(line, rect);
  }

}}}}  // namespace dials::algorithms::polygon::clip

#endif /* DIALS_ALGORITHMS_POLYGON_CLIPPING_CLIPPING_H */

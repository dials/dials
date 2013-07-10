/*
 * area.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_POLYGON_AREA_H
#define DIALS_ALGORITHMS_POLYGON_AREA_H

#include <scitbx/vec2.h>
#include <scitbx/array_family/flex_types.h>

namespace dials { namespace algorithms { namespace polygon {

  using scitbx::vec2;

  typedef scitbx::af::flex<vec2<double> >::type flex_vec2_double;

  /**
   * Calculate the area of a simple polygon using green's theorem
   * @param poly The polygon
   * @returns The area of the polygon
   */
  inline
  double simple_area(const flex_vec2_double &poly) {
    double a = 0.0;
    vec2<double> v0 = poly[poly.size()-1];
    for (std::size_t j = 0; j < poly.size(); ++j) {
      vec2<double> v1 = poly[j];
      a += v0[0] * v1[1] - v1[0] * v0[1];
      v0 = v1;
    }
    return a * 0.5;
  }

}}} // namespace dials::algorithms::polygon

#endif /* DIALS_ALGORITHMS_POLYGON_AREA_H */

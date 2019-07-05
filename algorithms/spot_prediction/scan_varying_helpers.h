/*
 * scan_varying_helpers.h
 *
 *  Copyright (C) 2013 CCP4, Diamond Light Source
 *
 *  Author: David Waterman (python code)
 *  Author: James Parkhurst (c++ port)
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */

#ifndef DIALS_ALGORITHMS_SPOT_PREDICTION_SCAN_VARYING_HELPERS_H
#define DIALS_ALGORITHMS_SPOT_PREDICTION_SCAN_VARYING_HELPERS_H

#include <cmath>
#include <scitbx/array_family/small.h>

namespace dials { namespace algorithms { namespace reeke_detail {

  /* Robust solution, for real roots only, of a quadratic in the form  (ax^2 +
   * bx + c)*/
  inline scitbx::af::small<double, 2> solve_quad(double a, double b, double c) {
    scitbx::af::small<double, 2> result;
    double discriminant = b * b - 4 * a * c;
    if (discriminant > 0.0) {
      int sign = (b >= 0 ? 1 : -1);
      double q = -0.5 * (b + sign * std::sqrt(discriminant));
      if (a != 0) result.push_back(q / a);
      if (q != 0) result.push_back(c / q);
    } else if (discriminant == 0.0) {
      double q = -b / (2 * a);
      result.push_back(q);
      result.push_back(q);
    }
    return result;
  }

}}}  // namespace dials::algorithms::reeke_detail

#endif  // DIALS_ALGORITHMS_SPOT_PREDICTION_SCAN_VARYING_HELPERS_H

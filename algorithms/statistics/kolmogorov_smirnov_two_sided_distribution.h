/*
 * kolmogorov_smirnov_two_sided_distribution.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */

#ifndef DIALS_ALGORITHMS_STATISTICS_KOLMOGOROV_SMIRNOV_TWO_SIDED_DISTRIBUTION_H
#define DIALS_ALGORITHMS_STATISTICS_KOLMOGOROV_SMIRNOV_TWO_SIDED_DISTRIBUTION_H

#include <iostream>
#include <limits>
#include <cmath>
#include <boost/math/special_functions.hpp>
#include <dials/error.h>

namespace dials { namespace algorithms {

  /**
   * The (asymptotic) kolmogorov_smirnov two sided distribution
   */
  template <typename RealType = double>
  class kolmogorov_smirnov_two_sided_distribution {
  public:
    typedef RealType value_type;
  };

  /**
   * An implementation of the limiting two sided kolmogorov smirnov distribution
   *
   * This implements the following equation.
   *
   *                  inf
   *  P(x) = 1 - 2 * SIGMA (-1)^(j-1) * exp(-2 * j^2 * x^2)
   *                  j=1
   *
   * @param dist The distribution
   * @param x A value between -inf and inf
   * @returns The value of the CDF at x
   */
  template <typename RealType>
  RealType cdf(const kolmogorov_smirnov_two_sided_distribution<RealType> &dist,
               const RealType &x) {
    const RealType epsilon = std::numeric_limits<RealType>::epsilon();
    RealType y = -2.0 * x * x;
    RealType s = 0.0, t = 0.0;
    int j = 1, sign = 1;
    do {
      t = std::exp(y * j * j);
      s += sign * t;
      ++j;
      sign *= -1;
    } while (s > 0.0 && (t / s) > epsilon);
    return 1.0 - 2.0 * s;
  }

}}  // namespace dials::algorithms

#endif  // DIALS_ALGORITHMS_STATISTICS_KOLMOGOROV_SMIRNOV_TWO_SIDED_DISTRIBUTION_H

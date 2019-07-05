/*
 * kolmogorov_smirnov_one_sided_distribution.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */

#ifndef DIALS_ALGORITHMS_STATISTICS_KOLMOGOROV_SMIRNOV_ONE_SIDED_DISTRIBUTION_H
#define DIALS_ALGORITHMS_STATISTICS_KOLMOGOROV_SMIRNOV_ONE_SIDED_DISTRIBUTION_H

#include <iostream>
#include <limits>
#include <cmath>
#include <boost/math/special_functions.hpp>
#include <dials/error.h>

namespace dials { namespace algorithms {

  /**
   * The kolmogorov smirnov one sided distribution
   */
  template <typename RealType = double>
  class kolmogorov_smirnov_one_sided_distribution {
  public:
    typedef RealType value_type;

    kolmogorov_smirnov_one_sided_distribution(std::size_t n) : n_(n) {
      DIALS_ASSERT(n > 0);
    }

    std::size_t n() const {
      return n_;
    }

  private:
    std::size_t n_;
  };

  namespace detail {

    /**
     * A direct implementation of the exact kolmogorov smirnov onesided
     * distribution as shown in
     *  Birnbaum, Z. & Tingey, F. (1951). Ann. Math. Stat. 592–596.
     *
     * This implements the following equation.
     *
     *                 n(1-x)
     *  P(x) = 1 - x * SIGMA (n) * (1-x-j/n)^(n-j) * (x - j/n)^(j-1)
     *                  j=0  (j)
     *
     * The direct implementation is reasonably fast but for large values of N,
     * it can result in a floating point overflow. It is valid for N where:
     *
     *  2^N < numeric_limits<RealType>::max
     */
    template <typename RealType>
    RealType cdf_small(const kolmogorov_smirnov_one_sided_distribution<RealType> &dist,
                       const RealType &x) {
      int n = (int)dist.n();
      int m = (int)std::floor(n * (1.0 - x));
      RealType s = 0.0;
      RealType b = 1.0;
      for (int j = 0; j <= m; ++j) {
        RealType a = x + (RealType)j / (RealType)n;
        s += b * std::pow(1.0 - a, n - j) * std::pow(a, j - 1);
        b *= (RealType)(n - j) / (RealType)(j + 1);
      }
      return 1.0 - x * s;
    }

    /**
     * An implementation fo the exact kolmogorov smirnov one sided distribution
     * (as shown above) for large values of N. This function is slower but uses
     * a sum of logarithms to avoid a floating point overflow. It uses the
     * relationship between the gamma function and the biniomial coefficients to
     * calculate log of the binomial coefficients
     */
    template <typename RealType>
    RealType cdf_large(const kolmogorov_smirnov_one_sided_distribution<RealType> &dist,
                       const RealType &x) {
      int n = (int)dist.n();
      int m = (int)std::floor(n * (1.0 - x));
      RealType s = 0.0;
      RealType b = boost::math::lgamma(n + 1);
      for (int j = 0; j <= m; ++j) {
        RealType a = x + (RealType)j / (RealType)n;
        if (1.0 - a > 0 && a > 0) {
          RealType c = b - boost::math::lgamma(j + 1) - boost::math::lgamma(n - j + 1)
                       + (n - j) * std::log(1.0 - a) + (j - 1) * std::log(a);
          s += std::exp(c);
        }
      }
      return 1.0 - x * s;
    }

    /**
     * A direct implementation of the PDF of the exact kolmogorov smirnov one
     * sided distribution as shown in Birnbaum, Z. & Tingey, F. (1951). Ann.
     * Math. Stat. 592–596.  This formula was dirived from the log form of the
     * sum.
     *
     * This implements the following equation.
     *
     *
     *             n(1-x)             n(1-x)
     *  dP(x) = - (SIGMA f(x,j) + x * SIGMA f(x,j) * (  (j - 1)  -  (n - j)))
     *              j=0                j=0             ln(x+j/n)   ln(1-x-j/n)
     *
     * where:
     *
     *  f(x, j) = (n) * (1-x-j/n)^(n-j) * (x - j/n)^(j-1)
     *            (j)
     *
     * The direct implementation is reasonably fast but for large values of N,
     * it can result in a floating point overflow. It is valid for N where:
     *
     *  2^N < numeric_limits<RealType>::max
     */
    // template <typename RealType>
    // RealType pdf_small(
    // const kolmogorov_smirnov_one_sided_distribution<RealType> &dist,
    // const RealType &x) {
    // int n = (int)dist.n();
    // int m = (int)std::floor(n * (1.0 - x));
    // RealType s1 = 0.0;
    // RealType s2 = 0.0;
    // RealType b = 1.0;
    // for (int j = 0; j <= m; ++j) {
    // RealType a = x + (RealType)j / (RealType)n;
    // RealType c = b * std::pow(1.0 - a, n - j) * std::pow(a, j - 1);
    // s1 += c;
    // RealType d = 1.0;
    ////        std::cout << a << std::endl;
    // if (a > 0 && a < 1.0) {
    // d = ((j - 1) / std::log(a) - (n - j) / std::log(1.0 - a));
    //} [>else if (a <= 0) {
    // d *= (n - j) / std::log(1.0);
    //} else if (a >= 1.0) {
    // d *= (j - 1) / std::log(1.0);
    //}*/
    // std::cout << x << ", " << d << std::endl;
    // s2 += c * d;
    ////if (a > 0 && a < 1.0) {
    ////s2 += c * ((j - 1) / std::log(a) - (n - j) / std::log(1.0 - a));
    ////}
    // b *= (RealType)(n - j) / (RealType)(j + 1);
    //}
    // return -(s1 + x * s2);
    //}

    /**
     * An implementation fo the exact kolmogorov smirnov one sided distribution
     * (as shown above) for large values of N. This function is slower but uses
     * a sum of logarithms to avoid a floating point overflow. It uses the
     * relationship between the gamma function and the biniomial coefficients to
     * calculate log of the binomial coefficients
     */
    // template <typename RealType>
    // RealType pdf_large(
    // const kolmogorov_smirnov_one_sided_distribution<RealType> &dist,
    // const RealType &x) {
    // int n = (int)dist.n();
    // int m = (int)std::floor(n * (1.0 - x));
    // RealType s1 = 0.0;
    // RealType s2 = 0.0;
    // RealType b = boost::math::lgamma(n + 1);
    // for (int j = 0; j <= m; ++j) {
    // RealType a = x + (RealType)j / (RealType)n;
    // RealType c = b
    //- boost::math::lgamma(j + 1) - boost::math::lgamma(n - j + 1)
    //+ (n - j)*std::log(1.0 - a) + (j - 1)*std::log(a);
    // RealType d = std::exp(c);
    // s1 += d;
    // s2 += d * ((j - 1) / std::log(a) - (n - j) / std::log(1.0 - a));
    //}
    // return -(s1 + x * s2);
    //}

  }  // namespace detail

  /**
   * Return the CDF for the kolmogorov smirnov one sided distribution.
   * @param dist The distribution
   * @param x A value between 0 and 1
   * @returns The value of the CDF at x
   */
  template <typename RealType>
  RealType cdf(const kolmogorov_smirnov_one_sided_distribution<RealType> &dist,
               const RealType &x) {
    DIALS_ASSERT(x >= 0 && x <= 1.0);
    DIALS_ASSERT(dist.n() > 0);
    if (x == 0.0) {
      return 0.0;
    } else if (x == 1.0) {
      return 1.0;
    }
    return dist.n() < std::numeric_limits<RealType>::max_exponent - 1
             ? detail::cdf_small(dist, x)
             : detail::cdf_large(dist, x);
  }

  /**
   * Return the PDF for the kolmogorov smirnov one sided distribution.
   * @param dist The distribution
   * @param x A value between 0 and 1
   * @returns The value of the PDF at x
   */
  // template <typename RealType>
  // RealType pdf(
  // const kolmogorov_smirnov_one_sided_distribution<RealType> &dist,
  // const RealType &x) {
  // DIALS_ASSERT(x >= 0 && x <= 1.0);
  // DIALS_ASSERT(dist.n() > 0);
  // if (x == 0.0) {
  // return 0.0;
  //} else if (x == 1.0) {
  // return 0.0;
  //}
  // return dist.n() < std::numeric_limits<RealType>::max_exponent - 1
  //? detail::pdf_small(dist, x)
  //: detail::pdf_large(dist, x);
  //}

}}  // namespace dials::algorithms

#endif  // DIALS_ALGORITHMS_STATISTICS_KOLMOGOROV_SMIRNOV_ONE_SIDED_DISTRIBUTION_H

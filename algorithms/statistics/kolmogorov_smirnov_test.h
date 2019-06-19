/*
 * kolmogorov_smirnov_test.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */

#ifndef DIALS_ALGORITHMS_STATISTICS_KOLMOGOROV_SMIRNOV_TEST_H
#define DIALS_ALGORITHMS_STATISTICS_KOLMOGOROV_SMIRNOV_TEST_H

#include <iostream>
#include <dials/algorithms/statistics/kolmogorov_smirnov_one_sided_distribution.h>
#include <dials/algorithms/statistics/kolmogorov_smirnov_two_sided_distribution.h>
#include <dials/error.h>

namespace dials { namespace algorithms {

  /**
   * Enumeration for type of test
   */
  enum KSType { Less, Greater, TwoSided };

  /**
   * Calculate D-
   */
  template <typename RealType>
  RealType kolmogorov_smirnov_test_d_minus(const std::vector<RealType> &cdfv) {
    RealType Dmax = 0.0;
    std::size_t n = cdfv.size();
    for (std::size_t i = 0; i < n; ++i) {
      RealType D = cdfv[i] - (RealType)i / (RealType)n;
      if (D > Dmax) {
        Dmax = D;
      }
    }
    return Dmax;
  }

  /**
   * Calculate D+
   */
  template <typename RealType>
  RealType kolmogorov_smirnov_test_d_plus(const std::vector<RealType> &cdfv) {
    RealType Dmax = 0.0;
    std::size_t n = cdfv.size();
    for (std::size_t i = 0; i < n; ++i) {
      RealType D = (RealType)(i + 1) / (RealType)n - cdfv[i];
      if (D > Dmax) {
        Dmax = D;
      }
    }
    return Dmax;
  }

  /**
   * Perform the test for one-sided less
   */
  template <typename RealType>
  std::pair<RealType, RealType> kolmogorov_smirnov_test_less(
    const std::vector<RealType> &cdfv) {
    typedef kolmogorov_smirnov_one_sided_distribution<RealType> ks_dist;
    RealType Dm = kolmogorov_smirnov_test_d_minus(cdfv);
    return std::make_pair(Dm, 1.0 - cdf(ks_dist(cdfv.size()), Dm));
  }

  /**
   * Perform the test for one-sided greater
   */
  template <typename RealType>
  std::pair<RealType, RealType> kolmogorov_smirnov_test_greater(
    const std::vector<RealType> &cdfv) {
    typedef kolmogorov_smirnov_one_sided_distribution<RealType> ks_dist;
    RealType Dp = kolmogorov_smirnov_test_d_plus(cdfv);
    return std::make_pair(Dp, 1.0 - cdf(ks_dist(cdfv.size()), Dp));
  }

  /**
   * Perform the two-sided test
   */
  template <typename RealType>
  std::pair<RealType, RealType> kolmogorov_smirnov_test_two_sided(
    const std::vector<RealType> &cdfv) {
    typedef kolmogorov_smirnov_one_sided_distribution<RealType> ks_dist1;
    typedef kolmogorov_smirnov_two_sided_distribution<RealType> ks_dist2;
    std::size_t n = cdfv.size();
    RealType Dm = kolmogorov_smirnov_test_d_minus(cdfv);
    RealType Dp = kolmogorov_smirnov_test_d_plus(cdfv);
    RealType D = std::max(Dm, Dp);
    RealType pa = 1.0 - cdf(ks_dist2(), (D * std::sqrt((RealType)n)));
    if (n > 2666 || pa > 0.8 - n * 0.3 / 1000.0) {
      return std::make_pair(D, pa);
    }
    return std::make_pair(D, (1.0 - cdf(ks_dist1(n), D)) * 2.0);
  }

  /**
   * Perform the kolmogorov smirnov test. Sorted the data into ascending order
   * the calculate the emiprical distribution function. Then compute the value
   * of the given distributions CDF at each sample value. The difference between
   * the empirical distribution and the CDF are then calculated and compared
   * against approximate values for the critical values of the kolmogorov
   * distribution.
   * @param dist The distribution
   * @param first The beginning of the data
   * @param last The end of the data
   * @param kstype The type of test to perform (less, greater, two_sided)
   * @returns (D, p-value)
   */
  template <typename Dist, typename Iterator>
  std::pair<typename Dist::value_type, typename Dist::value_type>
  kolmogorov_smirnov_test(const Dist &dist,
                          Iterator first,
                          Iterator last,
                          const KSType &kstype) {
    typedef typename Dist::value_type value_type;

    // Sort the sample values into ascending order
    std::vector<value_type> x(first, last);
    std::sort(x.begin(), x.end());

    // Calculate the value of the CDF at the sample points
    std::vector<value_type> cdfv(x.size());
    for (std::size_t i = 0; i < x.size(); ++i) {
      cdfv[i] = cdf(dist, x[i]);
    }

    // Do the ks test
    std::pair<value_type, value_type> result(0, 0);
    switch (kstype) {
    case Less:
      result = kolmogorov_smirnov_test_less(cdfv);
      break;
    case Greater:
      result = kolmogorov_smirnov_test_greater(cdfv);
      break;
    case TwoSided:
      result = kolmogorov_smirnov_test_two_sided(cdfv);
      break;
    default:
      DIALS_ASSERT(false);
      break;
    };
    return result;
  }

}}  // namespace dials::algorithms

#endif  // DIALS_ALGORITHMS_STATISTICS_KOLMOGOROV_SMIRNOV_TEST_H

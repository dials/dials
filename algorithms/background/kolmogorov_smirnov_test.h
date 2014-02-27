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
#ifndef DIALS_ALGORITHMS_BACKGROUND_KOLMOGOROV_SMIRNOV_TEST
#define DIALS_ALGORITHMS_BACKGROUND_KOLMOGOROV_SMIRNOV_TEST

#include <cmath>
#include <boost/math/special_functions/erf.hpp>
#include <scitbx/math/mean_and_variance.h>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/error.h>

namespace dials { namespace algorithms {

  using boost::math::erf;
  using scitbx::math::mean_and_variance;

  /**
   * Find the interval containing b
   * @param a The array of descending values
   * @param n The number of values
   * @param b The value to find
   * @returns The lower bound of the interval
   */
  inline
  std::size_t ks_critical_find_interval(double *a, std::size_t n, double b) {
    for (std::size_t i = 0; i < n - 1; ++i) {
      if (a[i] >= b && b >= a[i+1]) {
        return i;
      }
    }
    DIALS_ASSERT(false);
    return 0;
  }

  /**
   * Get the critical value to test against. These are approximate values.
   * @param n The number of points in the sample
   * @param a The probability (1 - a)
   * @returns The critical value
   */
  inline
  double ks_critical_large_sample(std::size_t n, double a) {

    // Setup the tables
    const std::size_t NALPHA = 5;
    double alpha[NALPHA] = { 0.20, 0.10, 0.05, 0.02, 0.01 };
    double apprx[NALPHA] = { 1.07, 1.22, 1.36, 1.52, 1.63 };

    // Find the lower bound of the interval
    std::size_t i1 = ks_critical_find_interval(alpha, NALPHA, a);

    // Compute the critical value. If N is less than 50 then
    // interpolate using the large sample approximation, otherwise,
    // lookup the values from the table and interpolate
    DIALS_ASSERT(n > 0);
    double sqrt_n = std::sqrt((double)n);
    double d1 = apprx[i1] / sqrt_n;
    double d2 = apprx[i1+1] / sqrt_n;
    double a1 = alpha[i1];
    double a2 = alpha[i1+1];
    return d1 + (d2 - d1) * (a - a1) / (a2 - a1);
  }

  /**
   * Get the critical value to test against. These are approximate values.
   * @param n The number of points in the sample
   * @param a The probability (1 - a)
   * @returns The critical value
   */
  inline
  double ks_critical_small_sample(std::size_t n, double a) {

    // Setup the tables
    const std::size_t NALPHA = 5;
    const std::size_t NAPPRX = 40;

    // The alpha values
    double alpha[NALPHA] = { 0.20, 0.10, 0.05, 0.02, 0.01 };

    // The apprixmate critical values
    double apprx[NAPPRX][NALPHA] = {
      { 0.900, 0.950, 0.975, 0.990, 0.995, },
      { 0.684, 0.776, 0.842, 0.900, 0.929, },
      { 0.565, 0.636, 0.708, 0.785, 0.829, },
      { 0.493, 0.565, 0.624, 0.689, 0.734, },
      { 0.447, 0.509, 0.563, 0.627, 0.669, },
      { 0.410, 0.468, 0.519, 0.577, 0.617, },
      { 0.381, 0.436, 0.483, 0.538, 0.576, },
      { 0.358, 0.410, 0.454, 0.507, 0.542, },
      { 0.339, 0.387, 0.430, 0.480, 0.513, },
      { 0.323, 0.369, 0.409, 0.457, 0.489, },
      { 0.308, 0.352, 0.391, 0.437, 0.468, },
      { 0.296, 0.338, 0.375, 0.419, 0.449, },
      { 0.285, 0.325, 0.361, 0.404, 0.432, },
      { 0.275, 0.314, 0.349, 0.390, 0.418, },
      { 0.266, 0.304, 0.338, 0.377, 0.404, },
      { 0.258, 0.295, 0.327, 0.366, 0.392, },
      { 0.250, 0.286, 0.318, 0.355, 0.381, },
      { 0.244, 0.279, 0.309, 0.346, 0.371, },
      { 0.237, 0.271, 0.301, 0.337, 0.361, },
      { 0.232, 0.265, 0.294, 0.329, 0.352, },
      { 0.226, 0.259, 0.287, 0.321, 0.344, },
      { 0.221, 0.253, 0.281, 0.314, 0.337, },
      { 0.216, 0.247, 0.275, 0.307, 0.330, },
      { 0.212, 0.242, 0.269, 0.301, 0.323, },
      { 0.208, 0.238, 0.264, 0.295, 0.317, },
      { 0.204, 0.233, 0.259, 0.290, 0.311, },
      { 0.200, 0.229, 0.254, 0.284, 0.305, },
      { 0.197, 0.225, 0.250, 0.279, 0.300, },
      { 0.193, 0.221, 0.246, 0.275, 0.295, },
      { 0.190, 0.218, 0.242, 0.270, 0.290, },
      { 0.187, 0.214, 0.238, 0.266, 0.285, },
      { 0.184, 0.211, 0.234, 0.262, 0.281, },
      { 0.182, 0.208, 0.231, 0.258, 0.277, },
      { 0.179, 0.205, 0.227, 0.254, 0.273, },
      { 0.177, 0.202, 0.224, 0.251, 0.269, },
      { 0.174, 0.199, 0.221, 0.247, 0.265, },
      { 0.172, 0.196, 0.218, 0.244, 0.262, },
      { 0.170, 0.194, 0.215, 0.241, 0.258, },
      { 0.168, 0.191, 0.213, 0.238, 0.255, },
      { 0.165, 0.189, 0.210, 0.235, 0.252, },
    };

    // Ensure n is valid
    DIALS_ASSERT(n > 0 && n <= 40);

    // Find the lower bound of the interval
    std::size_t i1 = ks_critical_find_interval(alpha, NALPHA, a);

    // Compute the critical value. If N is less than 50 then
    // interpolate using the large sample approximation, otherwise,
    // lookup the values from the table and interpolate
    double d1 = apprx[n-1][i1];
    double d2 = apprx[n-1][i1+1];
    double a1 = alpha[i1];
    double a2 = alpha[i1+1];
    return d1 + (d2 - d1) * (a - a1) / (a2 - a1);
  }

  /**
   * Get the critical value to test against. These are approximate values.
   * @param n The number of points in the sample
   * @param a The probability (1 - a)
   * @returns The critical value
   */
  inline
  double ks_critical(std::size_t n, double a) {
    return (n > 40 ?
        ks_critical_large_sample(n, a) :
        ks_critical_small_sample(n, a));
  }

  /**
   * Standardize the data so that we can compare it to a standard normal
   * distribution.
   * @param data The data to standardize
   * @returns The standardized data
   */
  inline
  af::shared<double> standardize(const af::const_ref<double> &data) {
    mean_and_variance<double> mv(data);
    double m = mv.mean();
    double s = mv.unweighted_sample_standard_deviation();
    af::shared<double> result(data.size());
    for (std::size_t i = 0; i < data.size(); ++i) {
      result[i] = (data[i] - m) / s;
    }
    return result;
  }

  /**
   * Using the kolmogorov smirnov test, check that the data is normally
   * distributed.
   * @param data The array of data
   * @param alpha The probability (1 - alpha)
   * @returns True/False The data is normally distributed
   */
  inline
  bool ks_is_normally_distributed(
      const af::const_ref<double> &data, double alpha) {

    // Standardize data (x - m) / sig
    af::shared<double> x = standardize(data);

    // Sort the data into acsending order
    std::sort(x.begin(), x.end());

    // Calculate the EDF
    af::shared<double> edf(x.size());
    for (std::size_t i = 0; i < x.size(); ++i) {
      edf[i] = (i + 1) / (double)x.size();
    }

    // Compute the maximum distance between each point and the normal cdf
    double dn = 0.0;
    double r2 = std::sqrt(2.0);
    for (std::size_t i = 0; i < data.size(); ++i) {
      double d = std::abs(edf[i] - (erf(x[i] / r2) + 1.0) / 2.0);
      if (d > dn) {
        dn = d;
      }
    }

    // Check the maximum distance is less than the critical value
    return dn < ks_critical(data.size(), alpha);
  }

}} // namespace dials::algorithms

#endif // DIALS_ALGORITHMS_BACKGROUND_KOLMOGOROV_SMIRNOV_TEST

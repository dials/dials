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
   * Get the critical value to test against. These are approximate values.
   * Tables from https://www.utdallas.edu/~herve/MolinAbdi1998-LillieforsTechReport.pdf
   * @param n The number of points in the sample
   * @param a The probability (1 - a)
   * @returns The critical value
   */
  inline
  double ks_critical(std::size_t n, double a) {

    const std::size_t NALPHA = 5;
    const std::size_t NSAMPLE = 39;

    // The alpha values
    double alpha[NALPHA] = { 0.01, 0.05, 0.1, 0.15, 0.2 };

    // Scaling values to use with large n approximation c / sqrt(n)
    double large_approx[NALPHA] = { 1.63, 1.36, 1.22, 1.14, 1.07 };

    // The values of n at which the values are specified
    std::size_t num[NSAMPLE] = {
       4,  5,  6,  7,  8,  9, 10, 11, 12, 13,
      14, 15, 16, 17, 18, 19, 20, 25, 30, 31,
      32, 33, 34, 35, 36, 37, 38, 39, 40, 41,
      42, 43, 44, 45, 46, 47, 48, 49, 50
    };

    // Table of small approximations
    double small_approx[NSAMPLE][NALPHA] = {
      { 0.3027, 0.3216, 0.3456, 0.3754, 0.4129 },
      { 0.2893, 0.3027, 0.3188, 0.3427, 0.3959 },
      { 0.2694, 0.2816, 0.2982, 0.3245, 0.3728 },
      { 0.2521, 0.2641, 0.2802, 0.3041, 0.3504 },
      { 0.2387, 0.2502, 0.2649, 0.2875, 0.3331 },
      { 0.2273, 0.2382, 0.2522, 0.2744, 0.3162 },
      { 0.2171, 0.2273, 0.2410, 0.2616, 0.3037 },
      { 0.2080, 0.2179, 0.2306, 0.2506, 0.2905 },
      { 0.2004, 0.2101, 0.2228, 0.2426, 0.2812 },
      { 0.1932, 0.2025, 0.2147, 0.2337, 0.2714 },
      { 0.1869, 0.1959, 0.2077, 0.2257, 0.2627 },
      { 0.1811, 0.1899, 0.2016, 0.2196, 0.2545 },
      { 0.1758, 0.1843, 0.1956, 0.2128, 0.2477 },
      { 0.1711, 0.1794, 0.1902, 0.2071, 0.2408 },
      { 0.1666, 0.1747, 0.1852, 0.2018, 0.2345 },
      { 0.1624, 0.1700, 0.1803, 0.1965, 0.2285 },
      { 0.1589, 0.1666, 0.1764, 0.1920, 0.2226 },
      { 0.1429, 0.1498, 0.1589, 0.1726, 0.2010 },
      { 0.1315, 0.1378, 0.1460, 0.1590, 0.1848 },
      { 0.1291, 0.1353, 0.1432, 0.1559, 0.1820 },
      { 0.1274, 0.1336, 0.1415, 0.1542, 0.1798 },
      { 0.1254, 0.1314, 0.1392, 0.1518, 0.1770 },
      { 0.1236, 0.1295, 0.1373, 0.1497, 0.1747 },
      { 0.1220, 0.1278, 0.1356, 0.1478, 0.1720 },
      { 0.1203, 0.1260, 0.1336, 0.1454, 0.1695 },
      { 0.1188, 0.1245, 0.1320, 0.1436, 0.1677 },
      { 0.1174, 0.1230, 0.1303, 0.1421, 0.1653 },
      { 0.1159, 0.1214, 0.1288, 0.1402, 0.1634 },
      { 0.1147, 0.1204, 0.1275, 0.1386, 0.1616 },
      { 0.1131, 0.1186, 0.1258, 0.1373, 0.1599 },
      { 0.1119, 0.1172, 0.1244, 0.1353, 0.1573 },
      { 0.1106, 0.1159, 0.1228, 0.1339, 0.1556 },
      { 0.1095, 0.1148, 0.1216, 0.1322, 0.1542 },
      { 0.1083, 0.1134, 0.1204, 0.1309, 0.1525 },
      { 0.1071, 0.1123, 0.1189, 0.1293, 0.1512 },
      { 0.1062, 0.1113, 0.1180, 0.1282, 0.1499 },
      { 0.1047, 0.1098, 0.1165, 0.1269, 0.1476 },
      { 0.1040, 0.1089, 0.1153, 0.1256, 0.1463 },
      { 0.1030, 0.1079, 0.1142, 0.1246, 0.1457 },
    };

    DIALS_ASSERT(n >= 4);
    DIALS_ASSERT(a >= 0.01 && a <= 0.2);

    // Find the 2 alpha values to interpolate within
    std::size_t i1 = NALPHA, i2 = NALPHA;
    for (std::size_t i = 1 ; i < NALPHA; ++i) {
      if (alpha[i] >= a) {
        i1 = i - 1;
        i2 = i;
      }
    }
    DIALS_ASSERT(i1 < NALPHA && i2 < NALPHA);

    // Compute the critical value. If N is less than 50 then
    // interpolate using the large sample approximation, otherwise,
    // lookup the values from the table and interpolate
    double result = 0.0;
    if (n > 50) {
      double sqrt_n = std::sqrt((double)n);
      double d1 = large_approx[i1] / sqrt_n;
      double d2 = large_approx[i2] / sqrt_n;
      double a1 = (double)alpha[i1];
      double a2 = (double)alpha[i2];
      result = d1 + (d2 - d1) * (a - a1) / (a2 - a1);
    } else {
      std::size_t j1 = 0, j2 = 1;
      for (std::size_t j = 1; j < NALPHA; ++j) {
        if (num[j] >= n) {
          j1 = j-1;
          j2 = j;
        }
      }
      double d11 = small_approx[j1][i1];
      double d12 = small_approx[j1][i2];
      double d21 = small_approx[j2][i1];
      double d22 = small_approx[j2][i2];
      double a1 = (double)alpha[i1];
      double a2 = (double)alpha[i2];
      double n1 = (double)num[j1];
      double n2 = (double)num[j2];
      result = (
          d11 * (a2 - a)*(n2 - n) +
          d21 * (a - a1)*(n2 - n) +
          d12 * (a2 - a)*(n - n1) +
          d22 * (a - a1)*(n - n1)
          ) / ((a2 - a1) * (n2 - n1));
    }
    return result;
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

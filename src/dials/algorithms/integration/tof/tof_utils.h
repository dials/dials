#ifndef DIALS_ALGORITHMS_INTEGRATION_TOF_UTILS_H
#define DIALS_ALGORITHMS_INTEGRATION_TOF_UTILS_H
#include <cmath>

namespace dials { namespace algorithms {

  inline double exp_safe(double x) {
    if (x > 700.0) return std::exp(700.0);
    if (x < -700.0) return std::exp(-700.0);
    return std::exp(x);
  }
  inline double erfc_safe(double x) {
    if (x > 10.0) return std::erfc(10.0);
    if (x < -10.0) return std::erfc(-10.0);
    return std::erfc(x);
  }

  inline double erfcx_safe(double x) {
    if (x > 25.0) {
      return 1.0 / (x * std::sqrt(M_PI));
    }
    return std::exp(x * x) * std::erfc(x);
  }

  inline bool is_finite_double(double x) {
    return std::isfinite(x);
  }

  inline double simpson_integrate(scitbx::af::const_ref<double> y,
                                  scitbx::af::const_ref<double> x) {
    /*
     * Composite 1/3 implementation
     * based on https://en.wikipedia.org/wiki/Simpson%27s_rule
     */

    const size_t n = y.size();
    // Min amount of data required
    if (n < 2) return 0.0;

    // Points must be odd for even number of intervals
    // If even use trapezoid on final point
    if (n % 2 == 0) {
      double res = 0.0;
      size_t m = n - 1;

      // Main loop
      for (size_t i = 0; i + 2 < m; i += 2) {
        double h0 = x[i + 1] - x[i];
        double h1 = x[i + 2] - x[i + 1];
        res += (x[i + 2] - x[i]) / 6.0 * (y[i] + 4.0 * y[i + 1] + y[i + 2]);
      }

      // Trapezoid for last interval (m-1,m)
      res += 0.5 * (x[n - 1] - x[n - 2]) * (y[n - 2] + y[n - 1]);
      return res;
    } else {
      // No need for Trapezoid on last point
      double res = 0.0;

      // Main loop
      for (size_t i = 0; i + 2 < n; i += 2) {
        res += (x[i + 2] - x[i]) / 6.0 * (y[i] + 4.0 * y[i + 1] + y[i + 2]);
      }
      return res;
    }
  }
}}  // namespace dials::algorithms
#endif /* DIALS_ALGORITHMS_INTEGRATION_TOF_UTILS */
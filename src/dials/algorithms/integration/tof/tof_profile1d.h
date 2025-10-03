#ifndef DIALS_ALGORITHMS_INTEGRATION_TOF_TOF_PROFILE1D_H
#define DIALS_ALGORITHMS_INTEGRATION_TOF_TOFPROFILE1D_H

#include <vector>
#include <algorithm>
#include <cmath>
#include <limits>
#include <iostream>
#include <array>
#include <cassert>
#include <dials/array_family/scitbx_shared_and_versa.h>

#include <eigen3/Eigen/Dense>
#include <eigen3/unsupported/Eigen/NonLinearOptimization>

/*
1D profile fitting method based on
Yano, N., Yamada, T., Hosoya, T. et al.
Application of profile fitting method to neutron time-of-flight protein
single crystal diffraction data collected at the iBIX. Sci Rep 6,
36628 (2016). https://doi.org/10.1038/srep36628
*/
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

  static scitbx::af::shared<double> profile1d_func(scitbx::af::const_ref<double> tof,
                                                   double A,
                                                   double alpha,
                                                   double beta,
                                                   double sigma,
                                                   double T_ph) {
    const size_t m = tof.size();
    scitbx::af::shared<double> out(m, 0.0);

    // (Numbers) refer to equations in https://doi.org/10.1038/srep36628

    double sigma2 = sigma * sigma;
    double sigma_sqrt = std::sqrt(2.0 * sigma2);
    double N = (alpha * beta) / (2.0 * (alpha + beta));  // (5)

    for (size_t i = 0; i < m; ++i) {
      double dT = tof[i] - T_ph;                             // (11)
      double u = alpha * 0.5 * (alpha * sigma2 + 2.0 * dT);  // (7)
      double v = beta * 0.5 * (beta * sigma2 - 2.0 * dT);    // (8)
      double y = (alpha * sigma2 + dT) / sigma_sqrt;         // (9)
      double z = (beta * sigma2 - dT) / sigma_sqrt;          // (10)

      // Stable evaluation with erfcx
      double term1 = std::exp(u - y * y) * erfcx_safe(y);
      double term2 = std::exp(v - z * z) * erfcx_safe(z);

      double val = A * N * (term1 + term2);  // (1)
      if (!std::isfinite(val)) val = 1e-12;
      out[i] = val;
    }
    return out;
  }

  // Eigen Functor wrapper
  template <typename _Scalar, int NX = Eigen::Dynamic, int NY = Eigen::Dynamic>
  struct Functor {
    typedef _Scalar Scalar;
    enum { InputsAtCompileTime = NX, ValuesAtCompileTime = NY };
    typedef Eigen::Matrix<Scalar, InputsAtCompileTime, 1> InputType;
    typedef Eigen::Matrix<Scalar, ValuesAtCompileTime, 1> ValueType;
    typedef Eigen::Matrix<Scalar, ValuesAtCompileTime, InputsAtCompileTime>
      JacobianType;

    int m_inputs, m_values;
    Functor() : m_inputs(InputsAtCompileTime), m_values(ValuesAtCompileTime) {}
    Functor(int inputs, int values) : m_inputs(inputs), m_values(values) {}
    int inputs() const {
      return m_inputs;
    }
    int values() const {
      return m_values;
    }
  };

  struct TOFProfileFunctor : Functor<double> {
    scitbx::af::const_ref<double> tof;
    scitbx::af::const_ref<double> y_norm;  // Assumed normalized
    std::array<double, 5> min_bounds;
    std::array<double, 5> max_bounds;

    TOFProfileFunctor(scitbx::af::const_ref<double> tof_,
                      scitbx::af::const_ref<double> y_norm_,
                      const std::array<double, 5>& minb,
                      const std::array<double, 5>& maxb)
        : Functor<double>(5, static_cast<int>(tof_.size())),
          tof(tof_),
          y_norm(y_norm_) {
      min_bounds = minb;
      max_bounds = maxb;
    }

    inline Eigen::VectorXd clamp_params(const Eigen::VectorXd& x) const {
      Eigen::VectorXd xc = x;
      for (int i = 0; i < x.size(); ++i) {
        xc[i] = std::min(std::max(x[i], min_bounds[i]), max_bounds[i]);
      }
      return xc;
    }

    int operator()(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) const {
      Eigen::VectorXd xc = clamp_params(x);
      double A = xc[0];
      double alpha = xc[1];
      double beta = xc[2];
      double sigma = xc[3];
      double T_ph = xc[4];

      scitbx::af::shared<double> model =
        profile1d_func(tof, A, alpha, beta, sigma, T_ph);
      assert(model.size() == static_cast<size_t>(values()));
      for (int i = 0; i < values(); ++i) {
        fvec[i] = y_norm[i] - model[i];
      }
      return 0;
    }

    // Numerical Jacobian (central differences)
    int df(const Eigen::VectorXd& x, Eigen::MatrixXd& J) const {
      const double eps = std::sqrt(std::numeric_limits<double>::epsilon());
      Eigen::VectorXd xc = clamp_params(x);
      int p = inputs();
      int m = values();

      Eigen::VectorXd f0(m);
      operator()(xc, f0);

      for (int j = 0; j < p; ++j) {
        Eigen::VectorXd xp = xc;
        Eigen::VectorXd xm = xc;
        double step = (std::abs(xc[j]) + 1.0) * eps;
        xp[j] += step;
        xm[j] -= step;
        Eigen::VectorXd fp(m), fm(m);
        operator()(xp, fp);
        operator()(xm, fm);
        for (int i = 0; i < m; ++i) {
          J(i, j) = (fp[i] - fm[i]) / (2.0 * step);
        }
      }
      return 0;
    }
  };

  class TOFProfile1D {
  public:
    scitbx::af::const_ref<double> tof;
    scitbx::af::const_ref<double> intensities;  // raw intensities
    scitbx::af::shared<double> y_norm;          // normalized intensities
    double intensity_max;

    // params
    double A, alpha, beta, sigma, T_ph;
    std::array<double, 5> min_bounds;
    std::array<double, 5> max_bounds;

    TOFProfile1D(scitbx::af::const_ref<double> tof_,
                 scitbx::af::const_ref<double> intensities_,
                 double A_,
                 double alpha_,
                 double beta_,
                 double T_ph_,
                 const std::array<double, 2> alpha_bounds,
                 const std::array<double, 2> beta_bounds)
        : tof(tof_),
          intensities(intensities_),
          A(A_),
          alpha(alpha_),
          beta(beta_),
          sigma(1.0),
          T_ph(T_ph_) {
      DIALS_ASSERT(tof.size() > 0);
      DIALS_ASSERT(tof.size() == intensities.size());

      // Ensure no negative values
      intensity_max = 1.0;
      if (!(intensities.size() == 0)) {
        intensity_max = *std::max_element(intensities.begin(), intensities.end());
        if (intensity_max <= 0.0) intensity_max = 1.0;
      }

      // build normalized y vector
      const size_t n = intensities.size();
      y_norm.resize(n);
      for (size_t i = 0; i < n; ++i) {
        double v = intensities[i];
        if (!is_finite_double(v)) v = 0.0;
        if (v < 0) v = 0.0;
        y_norm[i] = v / intensity_max;
      }

      // Set initial sigma from estimating peak width
      sigma = estimate_sigma_from_fwhm(tof, y_norm.const_ref());

      // Param bounds (A, alpha, beta, sigma, T_ph)
      min_bounds = {1.0, alpha_bounds[0], beta_bounds[0], sigma / 4.0, tof.front()};
      max_bounds = {1000.0 * intensity_max,
                    alpha_bounds[1],
                    beta_bounds[1],
                    std::max(100., sigma * 4.0),
                    tof.back()};

      // Sanity check
      DIALS_ASSERT(A >= min_bounds[0] && A <= max_bounds[0]);
      DIALS_ASSERT(alpha >= min_bounds[1] && alpha <= max_bounds[1]);
      DIALS_ASSERT(beta >= min_bounds[2] && beta <= max_bounds[2]);
      DIALS_ASSERT(sigma >= min_bounds[3] && sigma <= max_bounds[3]);
      DIALS_ASSERT(T_ph >= min_bounds[4] && T_ph <= max_bounds[4]);
    }

    scitbx::af::shared<double> result() const {
      scitbx::af::shared<double> m = profile1d_func(tof, A, alpha, beta, sigma, T_ph);
      for (auto& v : m)
        v *= intensity_max;
      return m;
    }

    double estimate_sigma_from_fwhm(scitbx::af::const_ref<double> tof,
                                    scitbx::af::const_ref<double> y) {
      /*
       * Estimate sigma using full width at half maximum of peak in y
       */

      // Not enough data
      DIALS_ASSERT(tof.size() >= 3);

      // locate peak
      size_t imax = std::distance(y.begin(), std::max_element(y.begin(), y.end()));
      double ymax = y[imax];

      // Negative peak
      DIALS_ASSERT(ymax > 0.0);

      double half_max = 0.5 * ymax;

      // Search left crossing
      double tL = tof.front();
      for (size_t i = imax; i-- > 0;) {
        if (y[i] <= half_max && y[i + 1] > half_max) {
          double t0 = tof[i], t1 = tof[i + 1];
          double y0 = y[i], y1 = y[i + 1];
          double frac = (half_max - y0) / (y1 - y0);
          tL = t0 + frac * (t1 - t0);
          break;
        }
      }

      // Search right crossing
      double tR = tof.back();
      for (size_t i = imax; i + 1 < y.size(); ++i) {
        if (y[i] > half_max && y[i + 1] <= half_max) {
          double t0 = tof[i], t1 = tof[i + 1];
          double y0 = y[i], y1 = y[i + 1];
          double frac = (half_max - y0) / (y1 - y0);
          tR = t0 + frac * (t1 - t0);
          break;
        }
      }

      // Full width at half maximum
      double fwhm = std::max(tR - tL, 0.0);
      DIALS_ASSERT(fwhm > 0.0);

      // 2.354520045 = approx(sqrt(2ln2))
      double sigma0 = fwhm / 2.354820045;

      // Unphysical sigma (large as at least one sample spacing)
      double mean_dt = (tof.back() - tof.front()) / std::max<size_t>(tof.size() - 1, 1);
      sigma0 = std::max(sigma0, mean_dt);
      return sigma0;
    }

    double calc_intensity() const {
      scitbx::af::shared<double> r = result();
      return simpson_integrate(r.const_ref(), tof);
    }

    double calc_variance(scitbx::af::const_ref<double> projected_variance) const {
      scitbx::af::shared<double> profile =
        profile1d_func(tof, A, alpha, beta, sigma, T_ph);
      DIALS_ASSERT(profile.size() == projected_variance.size());

      // Normalized profile divided by amplitude (A)
      scitbx::af::shared<double> p(profile.size());
      DIALS_ASSERT(A > 0.0);
      for (size_t i = 0; i < p.size(); ++i) {
        p[i] = profile[i] / A;
      }

      // Variance of amplitude A
      double denom = 0.0;
      for (size_t i = 0; i < p.size(); ++i) {
        double var = projected_variance[i];
        if (var <= 1e-7) continue;
        denom += (p[i] * p[i]) / var;
      }
      if (denom <= 0.0) return 0.0;
      double varA = 1.0 / denom;

      // Integral of profile shape
      double integral_p = simpson_integrate(p.const_ref(), tof);

      // Propagate to integrated intensity
      return varA * (integral_p * integral_p);
    }

    bool fit(int maxfev = 2000, double xtol = 1e-8, double ftol = 1e-8) {
      /*
       * Does least-squares minimization and updates A, alpha, beta, sigma, T_ph
       * Returns success of fitting
       */

      const int ndata = static_cast<int>(tof.size());

      if (ndata < 5) {
        // Not enough data
        return false;
      }

      TOFProfileFunctor functor(tof, y_norm.const_ref(), min_bounds, max_bounds);

      typedef Eigen::LevenbergMarquardt<TOFProfileFunctor, double> LM;
      LM lm(functor);

      // Initial parameter vector
      Eigen::VectorXd x(5);
      x[0] = A;
      x[1] = alpha;
      x[2] = beta;
      x[3] = sigma;
      x[4] = T_ph;

      lm.parameters.maxfev = maxfev;
      lm.parameters.xtol = xtol;
      lm.parameters.ftol = ftol;

      int result = lm.minimize(x);

      // Clamp final params to bounds
      for (int i = 0; i < 5; ++i) {
        x[i] = std::min(std::max(x[i], min_bounds[i]), max_bounds[i]);
      }

      // Update params
      A = x[0];
      alpha = x[1];
      beta = x[2];
      sigma = x[3];
      T_ph = x[4];

      // check convergence
      if (result == Eigen::LevenbergMarquardtSpace::ImproperInputParameters) {
        std::cerr << "LM improper input parameters\n";
        return false;
      }
      if (result < 0) {
        std::cerr << "LM failed with code " << result << "\n";
      }

      return true;
    }

    bool trust_result(double I_prf, double var_prf, double I_sum, double var_sum) {
      /*
       * Check for reasonable variance and if I/sigma is not too far from summation
       */

      if (var_prf < 1e-7) {
        return false;
      }

      if (var_sum < 1e-7) {
        return false;
      }

      double i_sigma_sum = I_sum / std::sqrt(var_sum);
      double i_sigma_prf = I_prf / std::sqrt(var_prf);
      double z_score = (i_sigma_sum - i_sigma_prf) / std::sqrt(var_sum + var_prf);
      return std::abs(z_score) < 0.25;
    }
  };

}}  // namespace dials::algorithms
#endif /* DIALS_ALGORITHMS_INTEGRATION_TOF_TOF_PROFILE1D */
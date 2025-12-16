#ifndef DIALS_ALGORITHMS_INTEGRATION_TOF_TOF_PROFILE1D_H
#define DIALS_ALGORITHMS_INTEGRATION_TOF_TOF_PROFILE1D_H

#include <algorithm>
#include <cmath>
#include <limits>
#include <array>
#include <cassert>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <random>
#include "tof_utils.h"

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

  /*
   * Holds params required for profile1d
   */
  struct TOFProfile1DParams {
    double A;
    double alpha;
    double alpha_min;
    double alpha_max;
    double beta;
    double beta_min;
    double beta_max;
    int n_restarts;              // number of attempts when fitting
    bool optimize_profile;       // If false the profile is generated with input params
    bool show_profile_failures;  // Prints debugging information

    TOFProfile1DParams(double A,
                       double alpha,
                       double alpha_min,
                       double alpha_max,
                       double beta,
                       double beta_min,
                       double beta_max,
                       int n_restarts,
                       bool optimize_profile,
                       bool show_profile_failures)

        : A(A),
          alpha(alpha),
          alpha_min(alpha_min),
          alpha_max(alpha_max),
          beta(beta),
          beta_min(beta_min),
          beta_max(beta_max),
          n_restarts(n_restarts),
          optimize_profile(optimize_profile),
          show_profile_failures(show_profile_failures) {}
  };

  static scitbx::af::shared<double> profile1d_func(scitbx::af::const_ref<double> tof,
                                                   double A,
                                                   double alpha,
                                                   double beta,
                                                   double sigma,
                                                   double T_ph) {
    /*
     * func used to generate the actual profile
     * (Numbers) refer to equations in https://doi.org/10.1038/srep36628
     */

    const size_t m = tof.size();
    scitbx::af::shared<double> out(m, 0.0);

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

  struct TOFProfileFunctor {
    scitbx::af::const_ref<double> tof;
    scitbx::af::const_ref<double> y_norm;  // Assumed normalized
    std::array<double, 5> min_bounds;      // parameter bounds
    std::array<double, 5> max_bounds;      // parameter bounds
    int num_data_points, num_params;

    TOFProfileFunctor(scitbx::af::const_ref<double> tof_,
                      scitbx::af::const_ref<double> y_norm_,
                      const std::array<double, 5>& minb,
                      const std::array<double, 5>& maxb)
        : tof(tof_), y_norm(y_norm_) {
      min_bounds = minb;
      max_bounds = maxb;
      num_data_points = tof.size();
      num_params = 5;
    }

    int values() const {
      return num_data_points;
    }

    int inputs() const {
      return num_params;
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
      assert(model.size() == num_data_points);
      for (int i = 0; i < num_data_points; ++i) {
        fvec[i] = y_norm[i] - model[i];
      }
      return 0;
    }

    int df(const Eigen::VectorXd& x, Eigen::MatrixXd& J) const {
      const double eps = 1e-5;
      Eigen::VectorXd xc = clamp_params(x);
      J.resize(num_data_points, num_params);

      for (int j = 0; j < num_params; ++j) {
        // Perturb param
        double delta = eps * std::max(1.0, std::abs(xc[j]));
        Eigen::VectorXd xp = xc, xm = xc;
        xp[j] += delta;
        xm[j] -= delta;

        Eigen::VectorXd xpc = clamp_params(xp);
        Eigen::VectorXd xmc = clamp_params(xm);

        double step = xpc[j] - xmc[j];
        if (std::abs(step) < 1e-14) {
          J.col(j).setZero();
          continue;
        }

        Eigen::VectorXd fp(num_data_points), fm(num_data_points);
        operator()(xpc, fp);
        operator()(xmc, fm);

        // Central difference
        J.col(j) = (fp - fm) / step;
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
    int n_restarts;

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
                 const std::array<double, 2> beta_bounds,
                 int n_restarts_)
        : tof(tof_),
          intensities(intensities_),
          A(A_),
          alpha(alpha_),
          beta(beta_),
          sigma(1.0),
          T_ph(T_ph_),
          n_restarts(n_restarts_) {
      DIALS_ASSERT(tof.size() > 0);
      DIALS_ASSERT(tof.size() == intensities.size());

      // Get max intensity
      intensity_max = 1.0;
      if (!(intensities.size() == 0)) {
        intensity_max = *std::max_element(intensities.begin(), intensities.end());
        if (intensity_max <= 0.0) intensity_max = 1.0;
      }

      // Get normalized y vector
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
      min_bounds = {1., alpha_bounds[0], beta_bounds[0], sigma / 4.0, tof.front()};

      max_bounds = {1e4 * intensity_max,
                    alpha_bounds[1],
                    beta_bounds[1],
                    std::max(100., sigma * 4.0),
                    tof.back()};

      // Sanity check params
      DIALS_ASSERT(A >= min_bounds[0] && A <= max_bounds[0]);
      DIALS_ASSERT(alpha >= min_bounds[1] && alpha <= max_bounds[1]);
      DIALS_ASSERT(beta >= min_bounds[2] && beta <= max_bounds[2]);
      DIALS_ASSERT(sigma >= min_bounds[3] && sigma <= max_bounds[3]);
      DIALS_ASSERT(T_ph >= min_bounds[4] && T_ph <= max_bounds[4]);
    }

    scitbx::af::shared<double> result() const {
      /*
       * Generates (unnormalized) profile for all positions in tof
       */

      scitbx::af::shared<double> m = profile1d_func(tof, A, alpha, beta, sigma, T_ph);
      for (auto& v : m)
        v *= intensity_max;
      return m;
    }

    double estimate_sigma_from_fwhm(scitbx::af::const_ref<double> tof,
                                    scitbx::af::const_ref<double> y) {
      /*
       * Estimates sigma param using full width at half maximum of peak in y
       */

      // Not enough data
      if (tof.size() >= 3) {
        return 1.0;
      }

      // locate peak
      size_t imax = std::distance(y.begin(), std::max_element(y.begin(), y.end()));
      double ymax = y[imax];

      // Negative peak
      if (ymax <= 0.0) {
        return 1.0;
      }

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

      // Unphysical sigma (check large as at least one sample spacing)
      double mean_dt = (tof.back() - tof.front()) / std::max<size_t>(tof.size() - 1, 1);
      sigma0 = std::max(sigma0, mean_dt);
      return sigma0;
    }

    double calc_intensity() const {
      /**
       * Get overall intensity with Simpsons rule then divide by mean_dt to
       * approximate summation scale
       */

      scitbx::af::shared<double> r = result();
      double mean_dt = (tof[tof.size() - 1] - tof[0]) / (tof.size() - 1);
      return simpson_integrate(r.const_ref(), tof) / mean_dt;
    }

    std::size_t get_max_profile_index() {
      /*
       * Returns the index of the max of the profile
       */

      auto profile_result = this->result();
      auto max_profile_it =
        std::max_element(profile_result.begin(), profile_result.end());
      std::size_t max_profile_index =
        std::distance(profile_result.begin(), max_profile_it);
      return max_profile_index;
    }

    bool fit(std::size_t max_sum_index,  // Peak index of the projected intensity
             bool show_profile_failures,
             int maxfev = 200,
             double xtol = 1e-8,
             double ftol = 1e-8) {
      /*
       * Least-squares minimization
       * Updates A, alpha, beta, sigma, T_ph
       * If fitting fails, params are perturbed n_restarts to find a solution
       */

      // Check enough data for fitting
      const int ndata = static_cast<int>(tof.size());
      if (ndata < 5) return false;

      TOFProfileFunctor functor(tof, y_norm.const_ref(), min_bounds, max_bounds);
      typedef Eigen::LevenbergMarquardt<TOFProfileFunctor, double> LM;

      auto run_single_fit = [&](const Eigen::VectorXd& x_init,
                                double& final_error) -> bool {
        LM lm(functor);
        lm.parameters.maxfev = maxfev;
        lm.parameters.xtol = xtol;
        lm.parameters.ftol = ftol;

        Eigen::VectorXd x = x_init;
        int result = lm.minimize(x);
        if (result < 0) return false;

        x = functor.clamp_params(x);

        // Compute residual norm
        Eigen::VectorXd fvec(functor.num_data_points);
        final_error = fvec.squaredNorm();

        // Update fitted parameters
        A = x[0];
        alpha = x[1];
        beta = x[2];
        sigma = x[3];
        T_ph = x[4];

        return true;
      };

      // First fit attempt
      Eigen::VectorXd x0(5);
      x0 << A, alpha, beta, sigma, T_ph;
      double fit_resid = std::numeric_limits<double>::infinity();
      bool success = run_single_fit(x0, fit_resid);
      std::size_t max_profile_index;
      double I_prf, I_var;

      if (success) {
        I_prf = this->calc_intensity();
        if (this->trust_result(fit_resid,
                               I_prf,
                               max_sum_index,
                               max_profile_index,
                               show_profile_failures)) {
          return true;
        }
      }

      // Initial fit failed, perturb initial params
      std::mt19937 rng(std::random_device{}());
      std::uniform_real_distribution<double> unit_dist(0.0, 1.0);

      for (int i = 0; i < n_restarts; ++i) {
        Eigen::VectorXd x_try(5);
        for (int j = 0; j < 5; ++j) {
          double span = max_bounds[j] - min_bounds[j];
          double rand_frac = (unit_dist(rng) - 0.5) * 0.4;
          double perturbed = x0[j] + rand_frac * span;
          x_try[j] = std::max(min_bounds[j], std::min(perturbed, max_bounds[j]));
        }

        // Attempt fit
        success = run_single_fit(x_try, fit_resid);
        if (!success) continue;

        I_prf = this->calc_intensity();
        max_profile_index = this->get_max_profile_index();
        if (this->trust_result(fit_resid,
                               I_prf,
                               max_sum_index,
                               max_profile_index,
                               show_profile_failures)) {
          return true;
        }
      }

      return false;
    }

    bool trust_result(double error,
                      double I_prf,
                      std::size_t max_sum_index,
                      std::size_t max_profile_index,
                      bool show_error = false) {
      /*
       * Tests to check the fit is reasonable
       */
      if (!std::isfinite(error) || error <= 0.0) {
        if (show_error) {
          std::cerr << "profile1d fitting failure: invalid error value (error=" << error
                    << ")\n";
        }
        return false;
      }

      // Check reasonable intensity
      if (I_prf < 1e-7) {
        if (show_error) {
          std::cerr << "profile1d fitting failure: profile intensity too small (I_prf="
                    << I_prf << ")\n";
        }
        return false;
      }

      // Check peak position close to data peak
      int peak_delta =
        std::abs(static_cast<int>(max_sum_index) - static_cast<int>(max_profile_index));
      if (peak_delta > 3) {
        if (show_error) {
          std::cerr << "profile1d fitting failure: peak index mismatch (max_sum_index="
                    << max_sum_index << ", max_profile_index=" << max_profile_index
                    << ", delta=" << peak_delta << ")\n";
        }
        return false;
      }

      // Check peak isn't very flat
      auto m = result();
      double max_val = *std::max_element(m.begin(), m.end());
      double mean_val = std::accumulate(m.begin(), m.end(), 0.0) / m.size();
      double contrast = (max_val - mean_val) / (max_val + 1e-12);
      if (contrast < 0.1) {
        if (show_error) {
          std::cerr << "profile1d fitting failure: insufficient peak contrast "
                    << "(contrast=" << contrast << ", max_val=" << max_val
                    << ", mean_val=" << mean_val << ")\n";
        }
        return false;
      }

      // Check correlation with data
      double profile_peak = 0.0, data_peak = 0.0;
      double num = 0.0, denom_y = 0.0, denom_m = 0.0;

      for (std::size_t i = 0; i < tof.size(); ++i) {
        double y = y_norm[i];
        double p = m[i] / intensity_max;

        if (i == 0 || y > data_peak) {
          data_peak = y;
        }
        if (i == 0 || p > profile_peak) {
          profile_peak = p;
        }

        num += y * p;
        denom_y += y * y;
        denom_m += p * p;
      }

      double corr = num / std::sqrt(denom_y * denom_m + 1e-12);
      if (corr < 0.9) {
        if (show_error) {
          std::cerr << "profile1d fitting failure: low correlation (corr=" << corr
                    << ")\n";
        }
        return false;
      }

      // Check peak height is within 10% of data peak
      double peak_diff = std::abs(profile_peak - data_peak);
      if (peak_diff > data_peak * 0.1) {
        if (show_error) {
          std::cerr << "profile1d fitting failure: peak height mismatch "
                    << "(profile_peak=" << profile_peak << ", data_peak=" << data_peak
                    << ", diff=" << peak_diff << ")\n";
        }
        return false;
      }

      return true;
    }
  };

  bool fit_profile1d(
    scitbx::af::const_ref<double> projected_intensity,
    scitbx::af::const_ref<double> tof_z,
    TOFProfile1DParams& profile_params,
    double& I_prf_out,
    boost::optional<scitbx::af::shared<double>> line_profile_out = boost::none,
    bool update_params = false) {
    /**
     * Wrapper for fitting a given reflection
     * If line_profile_out is provided the profile is returned at every
     * position in tof_z
     */

    // Get T_ph (peak position)
    auto max_it =
      std::max_element(projected_intensity.begin(), projected_intensity.end());
    size_t max_index = std::distance(projected_intensity.begin(), max_it);
    double T_ph = tof_z[max_index];

    // Fit profile
    const std::array<double, 2> alpha_bounds = {profile_params.alpha_min,
                                                profile_params.alpha_max};
    const std::array<double, 2> beta_bounds = {profile_params.beta_min,
                                               profile_params.beta_max};

    TOFProfile1D profile(tof_z,
                         projected_intensity,
                         profile_params.A,
                         profile_params.alpha,
                         profile_params.beta,
                         T_ph,
                         alpha_bounds,
                         beta_bounds,
                         profile_params.n_restarts);

    bool profile_success = true;
    if (profile_params.optimize_profile) {
      profile_success = profile.fit(max_index, profile_params.show_profile_failures);
    }

    if (profile_success) {
      if (update_params) {
        profile_params.alpha = profile.alpha;
        profile_params.beta = profile.beta;
        profile_params.A = profile.A;
      }
      double I_prf = profile.calc_intensity();
      auto profile_result = profile.result();
      DIALS_ASSERT(projected_intensity.size() == profile_result.size());

      I_prf_out = I_prf;

      if (!line_profile_out) {
        return profile_success;
      }

      scitbx::af::shared<double> line_profile = *line_profile_out;
      DIALS_ASSERT(line_profile.size() == profile_result.size());
      for (std::size_t i = 0; i < profile_result.size(); ++i) {
        line_profile[i] = profile_result[i];
      }
      return profile_success;
    }
    return false;
  }

}}  // namespace dials::algorithms
#endif /* DIALS_ALGORITHMS_INTEGRATION_TOF_TOF_PROFILE1D */

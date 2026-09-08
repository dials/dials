#ifndef DIALS_ALGORITHMS_INTEGRATION_TOF_TOF_PROFILE_1D_IC_H
#define DIALS_ALGORITHMS_INTEGRATION_TOF_TOF_PROFILE_1D_IC_H

#include <algorithm>
#include <cmath>
#include <limits>
#include <array>
#include <cassert>
#include <numeric>
#include <random>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include "tof_utils.h"
#include "tof_profile_3d_ic.h"  // for ic_profile()

#include <eigen3/Eigen/Dense>
#include <eigen3/unsupported/Eigen/NonLinearOptimization>

/*
1D profile fitting method using the Ikeda-Carpenter function

based on
Ikeda, S.I Carpenter J.M, Wide-energy-range, high-resolution measurements of
neutron pulse shapes of polyethylene moderators, Nucl. Instrum. Methods Phys.
Res. A. 239(3), 536-544, (1985) https://doi.org/10.1016/0168-9002(85)90033-6,

and Mantid's IkedaCarpenterPV peak function:
https://docs.mantidproject.org/nightly/fitting/fitfunctions/IkedaCarpenterPV.html
*/
namespace dials { namespace algorithms {

  /*
   * Holds params required for profile_1d_ic
   */
  struct TOFProfile1DICParams {
    double A, A_min, A_max;  // fast-neutron decay constant
    double B, B_min, B_max;  // slow-neutron decay constant
    double R, R_min, R_max;  // fast/slow neutron ratio
    double HatWidth;         // half-width of rectangular convolution kernel (ToF bins)
    double KConv;            // Gaussian convolution kernel decay rate (per ToF bin^2)
    int n_restarts;          // number of attempts when fitting
    bool optimize_profile;   // If false the profile is generated with input params
    bool show_profile_failures;  // Prints debugging information

    TOFProfile1DICParams(double A,
                         double A_min,
                         double A_max,
                         double B,
                         double B_min,
                         double B_max,
                         double R,
                         double R_min,
                         double R_max,
                         double HatWidth,
                         double KConv,
                         int n_restarts,
                         bool optimize_profile,
                         bool show_profile_failures)

        : A(A),
          A_min(A_min),
          A_max(A_max),
          B(B),
          B_min(B_min),
          B_max(B_max),
          R(R),
          R_min(R_min),
          R_max(R_max),
          HatWidth(HatWidth),
          KConv(KConv),
          n_restarts(n_restarts),
          optimize_profile(optimize_profile),
          show_profile_failures(show_profile_failures) {}
  };

  static scitbx::af::shared<double> profile1d_ic_func(scitbx::af::const_ref<double> tof,
                                                      double Scale,
                                                      double A,
                                                      double B,
                                                      double R,
                                                      double T0,
                                                      double HatWidth,
                                                      double KConv) {
    /*
     * func used to generate the actual profile
     */

    scitbx::af::shared<double> profile = ic_profile(tof, A, B, R, T0, HatWidth, KConv);
    const std::size_t m = profile.size();
    scitbx::af::shared<double> out(m, 0.0);
    for (std::size_t i = 0; i < m; ++i) {
      double val = Scale * profile[i];
      if (!std::isfinite(val)) val = 1e-12;
      out[i] = val;
    }
    return out;
  }

  struct IC1DProfileFunctor {
    scitbx::af::const_ref<double> tof;
    scitbx::af::const_ref<double> y_norm;  // Assumed normalized
    std::array<double, 7> min_bounds;      // parameter bounds
    std::array<double, 7> max_bounds;      // parameter bounds
    int num_data_points, num_params;

    IC1DProfileFunctor(scitbx::af::const_ref<double> tof_,
                       scitbx::af::const_ref<double> y_norm_,
                       const std::array<double, 7>& minb,
                       const std::array<double, 7>& maxb)
        : tof(tof_), y_norm(y_norm_) {
      min_bounds = minb;
      max_bounds = maxb;
      num_data_points = tof.size();
      num_params = 7;
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
      double Scale = xc[0];
      double A = xc[1];
      double B = xc[2];
      double R = xc[3];
      double T0 = xc[4];
      double HatWidth = xc[5];
      double KConv = xc[6];

      scitbx::af::shared<double> model =
        profile1d_ic_func(tof, Scale, A, B, R, T0, HatWidth, KConv);
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

  class TOFProfile1DIC {
  public:
    scitbx::af::const_ref<double> tof;
    scitbx::af::const_ref<double> intensities;  // raw intensities
    scitbx::af::shared<double> y_norm;          // normalized intensities
    double intensity_max;
    int n_restarts;

    // params
    double Scale, A, B, R, T0, HatWidth, KConv;
    std::array<double, 7> min_bounds;
    std::array<double, 7> max_bounds;

    TOFProfile1DIC(scitbx::af::const_ref<double> tof_,
                   scitbx::af::const_ref<double> intensities_,
                   double A_,
                   double B_,
                   double R_,
                   double T_ph_,
                   double HatWidth_,
                   double KConv_,
                   const std::array<double, 2>& A_bounds,
                   const std::array<double, 2>& B_bounds,
                   const std::array<double, 2>& R_bounds,
                   int n_restarts_)
        : tof(tof_),
          intensities(intensities_),
          A(A_),
          B(B_),
          R(R_),
          HatWidth(HatWidth_),
          KConv(KConv_),
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
      const std::size_t n = intensities.size();
      y_norm.resize(n);
      for (std::size_t i = 0; i < n; ++i) {
        double v = intensities[i];
        if (!is_finite_double(v)) v = 0.0;
        if (v < 0) v = 0.0;
        y_norm[i] = v / intensity_max;
      }

      Scale = 1.0;

      // The raw Ikeda-Carpenter function peaks at T0 + ~2/A, so
      // centre it so the observed peak (T_ph) lines up with the model peak.
      const double tof_min = tof.front();
      const double tof_max = tof.back();
      const double tof_range = std::abs(tof_max - tof_min);
      const double T0_min = tof_min - 2.0 * tof_range;

      const double T0_max =
        std::max(T0_min, tof_max - 2.0 / std::max(A_bounds[0], 1e-8));
      T0 = std::min(std::max(T_ph_ - 2.0 / std::max(A, 1e-8), T0_min), T0_max);

      // Param bounds (Scale, A, B, R, T0, HatWidth, KConv)
      min_bounds = {1e-6, A_bounds[0], B_bounds[0], R_bounds[0], T0_min, 0.5, 0.02};
      max_bounds = {1e8, A_bounds[1], B_bounds[1], R_bounds[1], T0_max, 20.0, 5.0};

      // Sanity check params
      DIALS_ASSERT(A >= min_bounds[1] && A <= max_bounds[1]);
      DIALS_ASSERT(B >= min_bounds[2] && B <= max_bounds[2]);
      DIALS_ASSERT(R >= min_bounds[3] && R <= max_bounds[3]);
      HatWidth = std::min(std::max(HatWidth, min_bounds[5]), max_bounds[5]);
      KConv = std::min(std::max(KConv, min_bounds[6]), max_bounds[6]);
    }

    scitbx::af::shared<double> result() const {
      /*
       * Generates unnormalised profile for all positions in tof
       */

      scitbx::af::shared<double> m =
        profile1d_ic_func(tof, Scale, A, B, R, T0, HatWidth, KConv);
      for (auto& v : m)
        v *= intensity_max;
      return m;
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
       * Updates Scale, A, B, R, T0, HatWidth, KConv
       * If fitting fails, params are perturbed n_restarts to find a solution
       */

      // Check enough data for fitting
      const int ndata = static_cast<int>(tof.size());
      if (ndata < 5) return false;

      IC1DProfileFunctor functor(tof, y_norm.const_ref(), min_bounds, max_bounds);
      typedef Eigen::LevenbergMarquardt<IC1DProfileFunctor, double> LM;

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
        functor(x, fvec);
        final_error = fvec.squaredNorm();

        // Update fitted parameters
        Scale = x[0];
        A = x[1];
        B = x[2];
        R = x[3];
        T0 = x[4];
        HatWidth = x[5];
        KConv = x[6];

        return true;
      };

      // First fit attempt
      Eigen::VectorXd x0(functor.num_params);
      x0 << Scale, A, B, R, T0, HatWidth, KConv;
      double fit_resid = std::numeric_limits<double>::infinity();
      bool success = run_single_fit(x0, fit_resid);
      std::size_t max_profile_index;
      double I_prf;

      if (success) {
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

      // Initial fit failed, perturb initial params
      std::mt19937 rng(std::random_device{}());
      std::uniform_real_distribution<double> unit_dist(0.0, 1.0);

      for (int i = 0; i < n_restarts; ++i) {
        Eigen::VectorXd x_try(functor.num_params);
        for (int j = 0; j < functor.num_params; ++j) {
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
          std::cerr << "profile_1d_ic fitting failure: invalid error value (error="
                    << error << ")\n";
        }
        return false;
      }

      // Check reasonable intensity
      if (I_prf < 1e-7) {
        if (show_error) {
          std::cerr
            << "profile_1d_ic fitting failure: profile intensity too small (I_prf="
            << I_prf << ")\n";
        }
        return false;
      }

      // Check peak position close to data peak
      int peak_delta =
        std::abs(static_cast<int>(max_sum_index) - static_cast<int>(max_profile_index));
      if (peak_delta > 3) {
        if (show_error) {
          std::cerr
            << "profile_1d_ic fitting failure: peak index mismatch (max_sum_index="
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
          std::cerr << "profile_1d_ic fitting failure: insufficient peak contrast "
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
          std::cerr << "profile_1d_ic fitting failure: low correlation (corr=" << corr
                    << ")\n";
        }
        return false;
      }

      // Check peak height is within 10% of data peak
      double peak_diff = std::abs(profile_peak - data_peak);
      if (peak_diff > data_peak * 0.1) {
        if (show_error) {
          std::cerr << "profile_1d_ic fitting failure: peak height mismatch "
                    << "(profile_peak=" << profile_peak << ", data_peak=" << data_peak
                    << ", diff=" << peak_diff << ")\n";
        }
        return false;
      }

      return true;
    }
  };

  bool fit_profile_1d_ic(
    scitbx::af::const_ref<double> projected_intensity,
    scitbx::af::const_ref<double> tof_z,
    TOFProfile1DICParams& profile_params,
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
    std::size_t max_index = std::distance(projected_intensity.begin(), max_it);
    double T_ph = tof_z[max_index];

    // Fit profile
    const std::array<double, 2> A_bounds = {profile_params.A_min, profile_params.A_max};
    const std::array<double, 2> B_bounds = {profile_params.B_min, profile_params.B_max};
    const std::array<double, 2> R_bounds = {profile_params.R_min, profile_params.R_max};

    TOFProfile1DIC profile(tof_z,
                           projected_intensity,
                           profile_params.A,
                           profile_params.B,
                           profile_params.R,
                           T_ph,
                           profile_params.HatWidth,
                           profile_params.KConv,
                           A_bounds,
                           B_bounds,
                           R_bounds,
                           profile_params.n_restarts);

    bool profile_success = true;
    if (profile_params.optimize_profile) {
      profile_success = profile.fit(max_index, profile_params.show_profile_failures);
    }

    if (profile_success) {
      if (update_params) {
        profile_params.A = profile.A;
        profile_params.B = profile.B;
        profile_params.R = profile.R;
        profile_params.HatWidth = profile.HatWidth;
        profile_params.KConv = profile.KConv;
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
#endif /* DIALS_ALGORITHMS_INTEGRATION_TOF_TOF_PROFILE_1D_IC_H */

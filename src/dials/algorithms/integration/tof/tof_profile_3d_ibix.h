#ifndef DIALS_ALGORITHMS_INTEGRATION_TOF_PROFILE_3D_IBIX_H
#define DIALS_ALGORITHMS_INTEGRATION_TOF_PROFILE_3D_IBIX_H

#include <array>
#include <algorithm>
#include <cmath>
#include <limits>
#include <random>
#include <boost/optional.hpp>
#include "tof_profile_3d_bvg.h"
#include "tof_utils.h"
#include <scitbx/constants.h>
#include <eigen3/Eigen/Dense>
#include <eigen3/unsupported/Eigen/NonLinearOptimization>
#include <dials/array_family/scitbx_shared_and_versa.h>

/*
3D profile fitting method using a separable
iBIX ToF shape (Yano et al.) × Bivariate Gaussian (BVG) model (XY),

based on
Yano, N., Yamada, T., Hosoya, T. et al.
Application of profile fitting method to neutron time-of-flight protein
single crystal diffraction data collected at the iBIX. Sci Rep 6,
36628 (2016). https://doi.org/10.1038/srep36628

ToF shape uses the same approach as tof_profile_1d_ibix.h, paired
with the same separable BVG spatial model used by tof_profile_3d_ic.h

*/

namespace dials { namespace algorithms {

  /*
   * Holds params required for TOFProfile3DIBIX
   */

  struct TOFProfile3DIBIXParams {
    // iBIX ToF shape parameters
    double alpha, alpha_min, alpha_max;  // fast decay constant
    double beta, beta_min, beta_max;     // slow decay constant
    double sigma_min, sigma_max;         // Gaussian instrument-resolution width

    // BVG parameters
    double SigX_min, SigX_max;
    double SigY_min, SigY_max;
    double SigP, SigP_min, SigP_max;  // spatial correlation

    // Control
    int n_restarts;              // number of attempts when fitting
    bool optimize_profile;       // If false the profile is generated with input params
    double max_drift_factor;     // Scales the per-reflection dXdt/dYdt bounds
                                 // relative to the shoebox's spatial extent
                                 // divided by its ToF extent
    bool show_profile_failures;  // Prints debugging information

    TOFProfile3DIBIXParams(double alpha_,
                           double alpha_min_,
                           double alpha_max_,
                           double beta_,
                           double beta_min_,
                           double beta_max_,
                           double sigma_min_,
                           double sigma_max_,
                           double SigX_min_,
                           double SigX_max_,
                           double SigY_min_,
                           double SigY_max_,
                           double SigP_,
                           double SigP_min_,
                           double SigP_max_,
                           int n_restarts_,
                           bool optimize_profile_,
                           double max_drift_factor_,
                           bool show_profile_failures_)
        : alpha(alpha_),
          alpha_min(alpha_min_),
          alpha_max(alpha_max_),
          beta(beta_),
          beta_min(beta_min_),
          beta_max(beta_max_),
          sigma_min(sigma_min_),
          sigma_max(sigma_max_),
          SigX_min(SigX_min_),
          SigX_max(SigX_max_),
          SigY_min(SigY_min_),
          SigY_max(SigY_max_),
          SigP(SigP_),
          SigP_min(SigP_min_),
          SigP_max(SigP_max_),
          n_restarts(n_restarts_),
          optimize_profile(optimize_profile_),
          max_drift_factor(max_drift_factor_),
          show_profile_failures(show_profile_failures_) {}
  };

  // iBIX ToF function with unit amplitude for each tof[i].
  // (Numbers) refer to equations in https://doi.org/10.1038/srep36628
  static scitbx::af::shared<double> ibix_raw(scitbx::af::const_ref<double> tof,
                                             double alpha,
                                             double beta,
                                             double sigma,
                                             double T_ph) {
    const std::size_t m = tof.size();
    scitbx::af::shared<double> out(m, 0.0);

    const double sigma2 = sigma * sigma;
    const double sigma_sqrt = std::sqrt(2.0 * sigma2);
    const double N = (alpha * beta) / (2.0 * (alpha + beta));  // (5)

    for (std::size_t i = 0; i < m; ++i) {
      const double dT = tof[i] - T_ph;                             // (11)
      const double u = alpha * 0.5 * (alpha * sigma2 + 2.0 * dT);  // (7)
      const double v = beta * 0.5 * (beta * sigma2 - 2.0 * dT);    // (8)
      const double y = (alpha * sigma2 + dT) / sigma_sqrt;         // (9)
      const double z = (beta * sigma2 - dT) / sigma_sqrt;          // (10)

      // Stable evaluation with erfcx
      const double term1 = std::exp(u - y * y) * erfcx_safe(y);
      const double term2 = std::exp(v - z * z) * erfcx_safe(z);

      const double val = N * (term1 + term2);  // (1), with amplitude A folded
                                               // into the fitted Scale instead
      out[i] = (std::isfinite(val) && val > 0.0) ? val : 0.0;
    }
    return out;
  }

  struct Profile3DIBIXFunctor {
    const scitbx::af::versa<vec3<double>, af::c_grid<3>>
      coords;  // (x,y,t) centred at peak
    const scitbx::af::versa<double, af::c_grid<3>> intensities;
    const scitbx::af::versa<double, af::c_grid<3>> background_variances;
    scitbx::af::shared<double> tof_axis;  // 1D ToF, length = accessor()[2]
    std::array<double, 9> min_bounds, max_bounds;
    int num_data_points, num_params;
    mutable double cached_Scale;
    mutable Eigen::VectorXd last_params;
    mutable bool cache_valid;

    /*
     * Parameter mapping:
     *   x[0] = log(alpha)   fast decay
     *   x[1] = log(beta)    slow decay
     *   x[2] = log(sigma)   instrument-resolution width
     *   x[3] = T_ph         time onset, relative to centred peak
     *   x[4] = log(SigX)    BVG sigma x
     *   x[5] = log(SigY)    BVG sigma y
     *   x[6] = SigP         BVG correlation (-1,1), clamped
     *   x[7] = dXdt         BVG centre drift in x per unit ToF (peak walk)
     *   x[8] = dYdt         BVG centre drift in y per unit ToF (peak walk)
     */

    Profile3DIBIXFunctor(
      const scitbx::af::versa<vec3<double>, af::c_grid<3>> coords_,
      const scitbx::af::versa<double, af::c_grid<3>> intensities_,
      const scitbx::af::versa<double, af::c_grid<3>> background_variances_,
      const std::array<double, 9>& minb,
      const std::array<double, 9>& maxb)
        : coords(coords_),
          intensities(intensities_),
          background_variances(background_variances_),
          cached_Scale(1.0),
          cache_valid(false) {
      min_bounds = minb;
      max_bounds = maxb;
      num_data_points = coords.size();
      num_params = 9;
      last_params =
        Eigen::VectorXd::Constant(num_params, std::numeric_limits<double>::quiet_NaN());

      const int nz = coords.accessor()[2];
      tof_axis.resize(nz);
      for (int iz = 0; iz < nz; ++iz)
        tof_axis[iz] = coords(0, 0, iz)[2];
    }

    int values() const {
      return num_data_points;
    }
    int inputs() const {
      return num_params;
    }

    Eigen::VectorXd clamp_params(const Eigen::VectorXd& x) const {
      Eigen::VectorXd xc = x;
      for (int i = 0; i < num_params; ++i)
        xc[i] = std::min(std::max(x[i], min_bounds[i]), max_bounds[i]);
      return xc;
    }

    void extract_params(const Eigen::VectorXd& xc,
                        double& alpha,
                        double& beta,
                        double& sigma,
                        double& T_ph,
                        double& SigX,
                        double& SigY,
                        double& SigP,
                        double& dXdt,
                        double& dYdt) const {
      alpha = std::exp(xc[0]);
      beta = std::exp(xc[1]);
      sigma = std::exp(xc[2]);
      T_ph = xc[3];
      SigX = std::exp(xc[4]);
      SigY = std::exp(xc[5]);
      SigP = xc[6];
      dXdt = xc[7];
      dYdt = xc[8];
    }

    double calc_scale(const scitbx::af::shared<double>& tof_shape,
                      double SigX,
                      double SigY,
                      double SigP,
                      double dXdt,
                      double dYdt) const {
      // Analytical amplitude via weighted least squares

      const int nx = coords.accessor()[0];
      const int ny = coords.accessor()[1];
      const int nz = coords.accessor()[2];
      double num = 0.0, den = 0.0;
      for (int ix = 0; ix < nx; ++ix) {
        for (int iy = 0; iy < ny; ++iy) {
          const double x0 = coords(ix, iy, 0)[0];
          const double y0 = coords(ix, iy, 0)[1];
          for (int iz = 0; iz < nz; ++iz) {
            const double dt = coords(ix, iy, iz)[2];
            const double bvg =
              bvg_func(x0 - dXdt * dt, y0 - dYdt * dt, SigX, SigY, SigP);
            const double p = tof_shape[iz] * bvg;
            const double obs = intensities(ix, iy, iz);
            double var_b = background_variances(ix, iy, iz);
            if (!std::isfinite(var_b) || var_b < 0.0) var_b = 0.0;
            const double var = std::max(std::abs(obs) + var_b, 1.0);
            const double w = 1.0 / var;
            num += w * obs * p;
            den += w * p * p;
          }
        }
      }
      return (den > 0.0) ? (num / den) : 0.0;
    }

    int operator()(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) const {
      const Eigen::VectorXd xc = clamp_params(x);
      double alpha, beta, sigma, T_ph, SigX, SigY, SigP, dXdt, dYdt;
      extract_params(xc, alpha, beta, sigma, T_ph, SigX, SigY, SigP, dXdt, dYdt);

      const scitbx::af::shared<double> tof_shape =
        ibix_raw(tof_axis.const_ref(), alpha, beta, sigma, T_ph);

      const int nx = coords.accessor()[0];
      const int ny = coords.accessor()[1];
      const int nz = coords.accessor()[2];

      if (!cache_valid || (xc - last_params).cwiseAbs().maxCoeff() > 1e-6) {
        cached_Scale = calc_scale(tof_shape, SigX, SigY, SigP, dXdt, dYdt);
        last_params = xc;
        cache_valid = true;
      }
      const double Scale = cached_Scale;

      fvec.resize(num_data_points);
      int count = 0;
      for (int ix = 0; ix < nx; ++ix) {
        for (int iy = 0; iy < ny; ++iy) {
          const double x0 = coords(ix, iy, 0)[0];
          const double y0 = coords(ix, iy, 0)[1];
          for (int iz = 0; iz < nz; ++iz) {
            const double dt = coords(ix, iy, iz)[2];
            const double bvg =
              bvg_func(x0 - dXdt * dt, y0 - dYdt * dt, SigX, SigY, SigP);
            const double model = Scale * tof_shape[iz] * bvg;
            const double obs = intensities(ix, iy, iz);
            double var_b = background_variances(ix, iy, iz);
            if (!std::isfinite(var_b) || var_b < 0.0) var_b = 0.0;
            const double var = std::max(std::abs(obs) + var_b, 1.0);
            const double sigma_w = std::sqrt(var);
            const double diff = (obs - model) / sigma_w;
            fvec[count++] = std::isfinite(diff) ? diff : 1e6;
          }
        }
      }
      return 0;
    }

    // Finite-difference Jacobian (forward differences)
    int df(const Eigen::VectorXd& x, Eigen::MatrixXd& J) const {
      J.resize(num_data_points, num_params);
      const double eps = 1e-5;
      const Eigen::VectorXd xc = clamp_params(x);
      Eigen::VectorXd f0(num_data_points);
      operator()(xc, f0);
      for (int j = 0; j < num_params; ++j) {
        Eigen::VectorXd xp = xc;
        const double delta = eps * (1.0 + std::abs(xc[j]));
        xp[j] += delta;
        const Eigen::VectorXd xpc = clamp_params(xp);
        Eigen::VectorXd fp(num_data_points);
        operator()(xpc, fp);
        const double step = xpc[j] - xc[j];
        if (std::abs(step) < 1e-8)
          J.col(j).setZero();
        else
          J.col(j) = (fp - f0) / step;
      }
      return 0;
    }
  };

  class TOFProfile3DIBIX {
  public:
    const scitbx::af::versa<vec3<double>, af::c_grid<3>> coords;
    const scitbx::af::versa<double, af::c_grid<3>> intensities;
    const scitbx::af::versa<double, af::c_grid<3>> background_variances;
    scitbx::af::versa<double, af::c_grid<3>> y_norm;
    int n_restarts;
    Eigen::VectorXd params;
    std::array<double, 9> min_bounds, max_bounds;
    double tof_range;
    boost::optional<Profile3DIBIXFunctor> functor;

    TOFProfile3DIBIX(
      scitbx::af::const_ref<vec3<double>, af::c_grid<3>> coords_,
      const scitbx::af::versa<double, af::c_grid<3>> intensities_,
      const scitbx::af::versa<double, af::c_grid<3>> background_variances_,
      double alpha,
      double beta,
      double SigP,
      const std::array<double, 2>& alpha_bounds,
      const std::array<double, 2>& beta_bounds,
      const std::array<double, 2>& sigma_bounds,
      const std::array<double, 2>& SigX_bounds,
      const std::array<double, 2>& SigY_bounds,
      const std::array<double, 2>& SigP_bounds,
      int n_restarts_,
      double max_drift_factor)
        : coords(get_rel_coords_3d(coords_, intensities_)),
          intensities(intensities_),
          background_variances(background_variances_),
          n_restarts(n_restarts_) {
      const int n = intensities.size();
      y_norm.resize(intensities.accessor());
      for (int i = 0; i < n; ++i) {
        const double v = intensities[i];
        y_norm[i] = (std::isfinite(v) && v > 0.0) ? v : 0.0;
      }

      const int nz = coords.accessor()[2];
      const double tof_min = coords(0, 0, 0)[2];
      const double tof_max = coords(0, 0, nz - 1)[2];
      tof_range = std::abs(tof_max - tof_min);

      // sigma init: estimate from the FWHM of the ToF-projected intensity
      scitbx::af::shared<double> tof_axis(nz), proj(nz, 0.0);
      for (int iz = 0; iz < nz; ++iz)
        tof_axis[iz] = coords(0, 0, iz)[2];
      for (std::size_t ix = 0; ix < coords.accessor()[0]; ++ix)
        for (std::size_t iy = 0; iy < coords.accessor()[1]; ++iy)
          for (int iz = 0; iz < nz; ++iz)
            proj[iz] += y_norm(ix, iy, iz);
      const double sigma_init = std::min(
        std::max(estimate_sigma_from_fwhm(tof_axis.const_ref(), proj.const_ref()),
                 sigma_bounds[0]),
        sigma_bounds[1]);

      // T_ph init: the iBIX ToF function peaks near T_ph, so centering places
      // the data peak at t=0
      const double T_ph_init = 0.0;
      const double T_ph_min = tof_min;
      const double T_ph_max = tof_max;

      const BVGSpatialInit sp = estimate_bvg_spatial_init(coords,
                                                          y_norm,
                                                          SigP,
                                                          SigX_bounds,
                                                          SigY_bounds,
                                                          SigP_bounds,
                                                          tof_range,
                                                          max_drift_factor);

      // Bounds array: (log_alpha, log_beta, log_sigma, T_ph, log_SigX,
      // log_SigY, SigP, dXdt, dYdt)
      min_bounds = {std::log(alpha_bounds[0]),
                    std::log(beta_bounds[0]),
                    std::log(sigma_bounds[0]),
                    T_ph_min,
                    std::log(SigX_bounds[0]),
                    std::log(SigY_bounds[0]),
                    SigP_bounds[0],
                    -sp.dXdt_max,
                    -sp.dYdt_max};
      max_bounds = {std::log(alpha_bounds[1]),
                    std::log(beta_bounds[1]),
                    std::log(sigma_bounds[1]),
                    T_ph_max,
                    std::log(SigX_bounds[1]),
                    std::log(SigY_bounds[1]),
                    SigP_bounds[1],
                    sp.dXdt_max,
                    sp.dYdt_max};

      params.resize(9);
      params[0] = std::log(std::min(std::max(alpha, alpha_bounds[0]), alpha_bounds[1]));
      params[1] = std::log(std::min(std::max(beta, beta_bounds[0]), beta_bounds[1]));
      params[2] = std::log(sigma_init);
      params[3] = T_ph_init;
      params[4] = std::log(sp.SigX_init);
      params[5] = std::log(sp.SigY_init);
      params[6] = sp.SigP_init;
      params[7] = 0.0;  // dXdt: no drift assumed unless the fit finds evidence
      params[8] = 0.0;  // dYdt

      functor.emplace(coords, y_norm, background_variances, min_bounds, max_bounds);
    }

    static double estimate_sigma_from_fwhm(scitbx::af::const_ref<double> tof,
                                           scitbx::af::const_ref<double> y) {
      /*
       * Estimates sigma param using full width at half maximum of peak in y
       */

      if (tof.size() < 3) return 1.0;

      const std::size_t imax =
        std::distance(y.begin(), std::max_element(y.begin(), y.end()));
      const double ymax = y[imax];
      if (ymax <= 0.0) return 1.0;

      const double half_max = 0.5 * ymax;

      double tL = tof.front();
      for (std::size_t i = imax; i-- > 0;) {
        if (y[i] <= half_max && y[i + 1] > half_max) {
          const double t0 = tof[i], t1 = tof[i + 1];
          const double y0 = y[i], y1 = y[i + 1];
          const double frac = (half_max - y0) / (y1 - y0);
          tL = t0 + frac * (t1 - t0);
          break;
        }
      }

      double tR = tof.back();
      for (std::size_t i = imax; i + 1 < y.size(); ++i) {
        if (y[i] > half_max && y[i + 1] <= half_max) {
          const double t0 = tof[i], t1 = tof[i + 1];
          const double y0 = y[i], y1 = y[i + 1];
          const double frac = (half_max - y0) / (y1 - y0);
          tR = t0 + frac * (t1 - t0);
          break;
        }
      }

      const double fwhm = std::max(tR - tL, 0.0);
      if (fwhm <= 0.0) return 1.0;

      // 2.354820045 = approx(sqrt(8ln2))
      double sigma0 = fwhm / 2.354820045;
      const double mean_dt =
        (tof.back() - tof.front()) / std::max<std::size_t>(tof.size() - 1, 1);
      sigma0 = std::max(sigma0, mean_dt);
      return sigma0;
    }

    // Param getters
    double get_alpha() const {
      return std::exp(params[0]);
    }
    double get_beta() const {
      return std::exp(params[1]);
    }
    double get_sigma() const {
      return std::exp(params[2]);
    }
    double get_T_ph() const {
      return params[3];
    }
    double get_SigX() const {
      return std::exp(params[4]);
    }
    double get_SigY() const {
      return std::exp(params[5]);
    }
    double get_SigP() const {
      return params[6];
    }
    double get_dXdt() const {
      return params[7];
    }
    double get_dYdt() const {
      return params[8];
    }

    scitbx::af::versa<double, af::c_grid<3>> result() const {
      const double alpha = get_alpha(), beta = get_beta(), sigma = get_sigma(),
                   T_ph = get_T_ph();
      const double SigX = get_SigX(), SigY = get_SigY(), SigP = get_SigP();
      const double dXdt = get_dXdt(), dYdt = get_dYdt();
      const scitbx::af::shared<double> tof_shape =
        ibix_raw(functor->tof_axis.const_ref(), alpha, beta, sigma, T_ph);
      const int nx = coords.accessor()[0];
      const int ny = coords.accessor()[1];
      const int nz = coords.accessor()[2];
      const double Scale = functor->cached_Scale;
      scitbx::af::versa<double, af::c_grid<3>> out(coords.accessor());
      for (int ix = 0; ix < nx; ++ix) {
        for (int iy = 0; iy < ny; ++iy) {
          const double x0 = coords(ix, iy, 0)[0];
          const double y0 = coords(ix, iy, 0)[1];
          for (int iz = 0; iz < nz; ++iz) {
            const double dt = coords(ix, iy, iz)[2];
            const double bvg =
              bvg_func(x0 - dXdt * dt, y0 - dYdt * dt, SigX, SigY, SigP);
            const double val = Scale * tof_shape[iz] * bvg;
            out(ix, iy, iz) = std::isfinite(val) ? val : 0.0;
          }
        }
      }
      return out;
    }

    double calc_intensity() const {
      const scitbx::af::versa<double, af::c_grid<3>> r = result();
      double total = 0.0;
      for (std::size_t i = 0; i < r.size(); ++i)
        total += r[i];
      return total;
    }

    bool fit(bool show_profile_failures,
             int maxfev = 500,
             double xtol = 1e-8,
             double ftol = 1e-8) {
      typedef Eigen::LevenbergMarquardt<Profile3DIBIXFunctor, double> LM;

      auto run_single_fit = [&](const Eigen::VectorXd& x_init,
                                double& final_error) -> bool {
        Eigen::VectorXd x = x_init;
        LM lm(*functor);
        lm.parameters.maxfev = maxfev;
        lm.parameters.xtol = xtol;
        lm.parameters.ftol = ftol;
        const int result = lm.minimize(x);
        if (result < 0) return false;
        x = functor->clamp_params(x);
        Eigen::VectorXd fvec(functor->num_data_points);
        (*functor)(x, fvec);
        final_error = fvec.squaredNorm();
        params = x;
        return true;
      };

      const Eigen::VectorXd x0 = params;
      double fit_resid = std::numeric_limits<double>::infinity();
      bool success = run_single_fit(x0, fit_resid);

      if (success) {
        const double I_prf = this->calc_intensity();
        if (this->trust_result(fit_resid, I_prf, show_profile_failures)) return true;
      }

      std::mt19937 rng(std::random_device{}());
      std::uniform_real_distribution<double> unit_dist(0.0, 1.0);
      std::normal_distribution<double> norm_dist(0.0, 0.5);

      // Log-scale parameter indices: alpha(0), beta(1), sigma(2), SigX(4), SigY(5).
      const int n_log = 5;
      const int log_idx[5] = {0, 1, 2, 4, 5};

      // T_ph (position), SigP and dXdt/dYdt (peak walk) are not log-scale
      // params, so perturb them additively, scaled by their own bound width,
      // in every tier.
      auto perturb_position = [&](Eigen::VectorXd& x_try, double scale) {
        for (int j : {3, 6, 7, 8}) {
          x_try[j] += norm_dist(rng) * scale * (max_bounds[j] - min_bounds[j]);
        }
      };

      for (int i = 0; i < n_restarts; ++i) {
        Eigen::VectorXd x_try = x0;

        if (i < n_restarts / 3) {
          // Small log-scale perturbations
          for (int j = 0; j < n_log; ++j)
            x_try[log_idx[j]] += std::log(0.8 + 0.4 * unit_dist(rng));
          perturb_position(x_try, 0.1);
        } else if (i < 2 * n_restarts / 3) {
          // Larger uniform log scale
          const double scale = 0.5 + 1.5 * unit_dist(rng);
          for (int j = 0; j < n_log; ++j)
            x_try[log_idx[j]] += std::log(scale);
          perturb_position(x_try, 0.3);
        } else {
          // Random within bounds
          for (int j = 0; j < functor->num_params; ++j)
            x_try[j] = min_bounds[j] + unit_dist(rng) * (max_bounds[j] - min_bounds[j]);
        }
        x_try = functor->clamp_params(x_try);

        success = run_single_fit(x_try, fit_resid);
        if (!success) continue;
        const double I_prf = this->calc_intensity();
        if (this->trust_result(fit_resid, I_prf, show_profile_failures)) return true;
      }
      return false;
    }

    bool trust_result(double error, double I_prf, bool show_error = false) const {
      if (!std::isfinite(error) || error <= 0.0) {
        if (show_error)
          std::cerr << "profile_3d_ibix fitting failure: invalid error (" << error
                    << ")\n";
        return false;
      }
      if (I_prf < 1e-7) {
        if (show_error)
          std::cerr << "profile_3d_ibix fitting failure: intensity too small (" << I_prf
                    << ")\n";
        return false;
      }
      const double SigX = get_SigX(), SigY = get_SigY(), SigP = get_SigP();
      if (SigX <= 0.0 || SigY <= 0.0) {
        if (show_error)
          std::cerr << "profile_3d_ibix fitting failure: non-positive sigma\n";
        return false;
      }
      if (std::abs(SigP) >= 1.0) {
        if (show_error)
          std::cerr << "profile_3d_ibix fitting failure: |SigP|>=1 (" << SigP << ")\n";
        return false;
      }

      const scitbx::af::versa<double, af::c_grid<3>> r = result();
      const double limit = 1e12;
      double sum_d = 0.0, sum_m = 0.0;
      double profile_peak = 0.0, data_peak = 0.0;
      const std::size_t n_vox = r.size();

      for (std::size_t ix = 0; ix < coords.accessor()[0]; ++ix) {
        for (std::size_t iy = 0; iy < coords.accessor()[1]; ++iy) {
          for (std::size_t iz = 0; iz < coords.accessor()[2]; ++iz) {
            const double m = r(ix, iy, iz);
            if (!std::isfinite(m)) {
              if (show_error)
                std::cerr
                  << "profile_3d_ibix fitting failure: non-finite model value\n";
              return false;
            }
            if (m > limit) {
              if (show_error)
                std::cerr
                  << "profile_3d_ibix fitting failure: model value exceeds limit\n";
              return false;
            }
            const double d = y_norm(ix, iy, iz);
            data_peak = std::max(data_peak, d);
            profile_peak = std::max(profile_peak, m);
            sum_d += d;
            sum_m += m;
          }
        }
      }

      const double mean_d = sum_d / n_vox;
      const double mean_m = sum_m / n_vox;
      double num = 0.0, denom_d = 0.0, denom_m = 0.0;
      for (std::size_t ix = 0; ix < coords.accessor()[0]; ++ix) {
        for (std::size_t iy = 0; iy < coords.accessor()[1]; ++iy) {
          for (std::size_t iz = 0; iz < coords.accessor()[2]; ++iz) {
            const double d = y_norm(ix, iy, iz) - mean_d;
            const double m = r(ix, iy, iz) - mean_m;
            num += d * m;
            denom_d += d * d;
            denom_m += m * m;
          }
        }
      }

      const double corr = num / std::sqrt(denom_d * denom_m + 1e-12);
      if (corr < 0.75) {
        if (show_error)
          std::cerr << "profile_3d_ibix fitting failure: low correlation (" << corr
                    << ")\n";
        return false;
      }
      if (data_peak > 0.0 && std::abs(profile_peak - data_peak) > 0.25 * data_peak) {
        if (show_error)
          std::cerr << "profile_3d_ibix fitting failure: peak height mismatch (profile="
                    << profile_peak << ", data=" << data_peak << ")\n";
        return false;
      }
      return true;
    }
  };

  bool fit_profile_3d_ibix(
    scitbx::af::const_ref<vec3<double>, af::c_grid<3>> coords,
    const scitbx::af::versa<double, af::c_grid<3>> intensities,
    const scitbx::af::versa<double, af::c_grid<3>> background_variances,
    TOFProfile3DIBIXParams& profile_params,
    double& I_prf_out,
    boost::optional<scitbx::af::versa<double, af::c_grid<3>>> profile_3d_out =
      boost::none,
    bool update_params = false) {
    const std::array<double, 2> alpha_bounds = {profile_params.alpha_min,
                                                profile_params.alpha_max};
    const std::array<double, 2> beta_bounds = {profile_params.beta_min,
                                               profile_params.beta_max};
    const std::array<double, 2> sigma_bounds = {profile_params.sigma_min,
                                                profile_params.sigma_max};
    const std::array<double, 2> SigX_bounds = {profile_params.SigX_min,
                                               profile_params.SigX_max};
    const std::array<double, 2> SigY_bounds = {profile_params.SigY_min,
                                               profile_params.SigY_max};
    const std::array<double, 2> SigP_bounds = {profile_params.SigP_min,
                                               profile_params.SigP_max};

    TOFProfile3DIBIX profile(coords,
                             intensities,
                             background_variances,
                             profile_params.alpha,
                             profile_params.beta,
                             profile_params.SigP,
                             alpha_bounds,
                             beta_bounds,
                             sigma_bounds,
                             SigX_bounds,
                             SigY_bounds,
                             SigP_bounds,
                             profile_params.n_restarts,
                             profile_params.max_drift_factor);

    bool profile_success = true;
    if (profile_params.optimize_profile)
      profile_success = profile.fit(profile_params.show_profile_failures);

    if (profile_success) {
      I_prf_out = profile.calc_intensity();

      if (update_params) {
        profile_params.alpha = profile.get_alpha();
        profile_params.beta = profile.get_beta();
        profile_params.SigP = profile.get_SigP();
      }

      if (profile_3d_out) {
        const scitbx::af::versa<double, af::c_grid<3>> pred = profile.result();
        scitbx::af::versa<double, af::c_grid<3>>& out = *profile_3d_out;
        DIALS_ASSERT(pred.accessor().all_eq(out.accessor()));
        for (std::size_t x = 0; x < pred.accessor()[0]; ++x)
          for (std::size_t y = 0; y < pred.accessor()[1]; ++y)
            for (std::size_t z = 0; z < pred.accessor()[2]; ++z)
              out(x, y, z) = pred(x, y, z);
      }
    }
    return profile_success;
  }

}}  // namespace dials::algorithms

#endif  // DIALS_ALGORITHMS_INTEGRATION_TOF_PROFILE_3D_IBIX_H

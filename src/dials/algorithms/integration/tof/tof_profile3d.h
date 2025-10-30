#ifndef DIALS_ALGORITHMS_INTEGRATION_TOF_GUTMANNPROFILE3D_H
#define DIALS_ALGORITHMS_INTEGRATION_TOF_GUTMANNPROFILE3D_H

#include <vector>
#include <array>
#include <algorithm>
#include <cmath>
#include <limits>
#include <iostream>
#include <random>
#include <cassert>
#include "tof_utils.h"

#include <eigen3/Eigen/Dense>
#include <eigen3/unsupported/Eigen/NonLinearOptimization>
#include <dials/array_family/scitbx_shared_and_versa.h>

/*
3D profile fitting method based on
Gutmann, M.J, A 3D profile function suitable for integration of
neutron time-of-flight single crystal diffraction peaks,
Nucl. Instrum. Methods Phys. Res. A, 848, 170-173, (2016)
https://doi.org/10.1016/j.nima.2016.12.026
*/

namespace dials { namespace algorithms {

  /*
   * Holds params required for profile3d
   */

  struct TOFProfile3DParams {
    double alpha;
    double alpha_min;
    double alpha_max;
    double beta;
    double beta_min;
    double beta_max;
    int n_restarts;         // number of attempts when fitting
    bool optimize_profile;  // If true the profile is generated with input params

    TOFProfile3DParams(double alpha,
                       double alpha_min,
                       double alpha_max,
                       double beta,
                       double beta_min,
                       double beta_max,
                       int n_restarts,
                       bool optimize_profile)
        : alpha(alpha),
          alpha_min(alpha_min),
          alpha_max(alpha_max),
          beta(beta),
          beta_min(beta_min),
          beta_max(beta_max),
          n_restarts(n_restarts),
          optimize_profile(optimize_profile) {}
  };

  struct GutmannProfileFunctor {
    const scitbx::af::versa<vec3<double>, af::c_grid<3>> coords;
    scitbx::af::shared<double> dt_widths;
    const scitbx::af::versa<double, af::c_grid<3>> intensities;
    const scitbx::af::versa<double, af::c_grid<3>> background_variances;
    mutable Eigen::VectorXd last_params;  // size 3 or full param vector
    mutable double cached_norm = 1.0;
    mutable double cached_A = 1.0;
    mutable bool cache_valid = false;
    std::array<double, 8> min_bounds;
    std::array<double, 8> max_bounds;
    int m, n;

    GutmannProfileFunctor(
      const scitbx::af::versa<vec3<double>, af::c_grid<3>> coords_,
      const scitbx::af::versa<double, af::c_grid<3>> intensities_,
      const scitbx::af::versa<double, af::c_grid<3>> background_variances_,
      const std::array<double, 8>& minb,
      const std::array<double, 8>& maxb)
        : coords(coords_),
          intensities(intensities_),
          background_variances(background_variances_) {
      min_bounds = minb;
      max_bounds = maxb;
      m = coords.size();  // Num data points
      n = 8;              // Num params

      dt_widths.resize(coords_.accessor()[2]);
      dt_widths[0] = 1.;
      for (std::size_t c_z = 1; c_z < coords.accessor()[2]; ++c_z) {
        dt_widths[c_z] = coords(0, 0, c_z)[2] - coords(0, 0, c_z - 1)[2];
      }
    }

    int inputs() const {
      return n;
    }
    int values() const {
      return m;
    }

    inline Eigen::VectorXd clamp_params(const Eigen::VectorXd& x) const {
      Eigen::VectorXd xc = x;
      for (int i = 0; i < x.size(); ++i)
        xc[i] = std::min(std::max(x[i], min_bounds[i]), max_bounds[i]);
      return xc;
    }

    Eigen::Matrix3d build_H_from_L(const Eigen::VectorXd& params) const {
      double l11 = std::exp(params[0]);
      double l21 = params[1];
      double l31 = params[2];
      double l22 = std::exp(params[3]);
      double l32 = params[4];
      double l33 = std::exp(params[5]);

      Eigen::Matrix3d L;
      L << l11, 0.0, 0.0, l21, l22, 0.0, l31, l32, l33;

      return L * L.transpose();
    }

    double func(const vec3<double>& c,
                const Eigen::Matrix3d& H,
                double alpha,
                double beta,
                double A_,
                double norm_factor,
                double dt_width) const {
      double dx = c[0], dy = c[1], dt = c[2];

      double H1 = H(0, 0), H2 = H(0, 1), H3 = H(0, 2);
      double H4 = H(1, 1), H5 = H(1, 2), H6 = H(2, 2);

      double a = alpha;
      double b = beta;

      double N = (a * b) / (2.0 * (a + b));
      double detH = std::max(H.determinant(), 1e-12);
      double Ng = std::sqrt(detH) / std::pow(2.0 * M_PI, 1.5);
      double f1 = N * Ng * std::sqrt(M_PI / (2.0 * H6));

      double u = 0.5 * a * (a + 2.0 * H6 * dt + 2.0 * H3 * dx + 2.0 * H5 * dy);
      double v = 0.5 * b * (b - 2.0 * H6 * dt - 2.0 * H3 * dx - 2.0 * H5 * dy);

      v = std::max(-700.0, std::min(v, 700.0));
      u = std::max(-700.0, std::min(u, 700.0));

      double y = (a + H6 * dt + H3 * dx + H5 * dy) / std::sqrt(2.0 * H6);
      double w = (b - H6 * dt - H3 * dx - H5 * dy) / std::sqrt(2.0 * H6);

      double f2 =
        std::exp(-0.5 * H1 * dx * dx - H2 * dx * dy - 0.5 * H4 * dy * dy
                 + (H3 * H3 * dx * dx + 2.0 * H3 * H5 * dx * dy + H5 * H5 * dy * dy)
                     / (2.0 * H6));

      auto F_tof = [&](double dt) {
        double u =
          0.5 * alpha * (alpha + 2.0 * H6 * dt + 2.0 * H3 * dx + 2.0 * H5 * dy);
        double v = 0.5 * beta * (beta - 2.0 * H6 * dt - 2.0 * H3 * dx - 2.0 * H5 * dy);

        double y = (alpha + H6 * dt + H3 * dx + H5 * dy) / std::sqrt(2.0 * H6);
        double w = (beta - H6 * dt - H3 * dx - H5 * dy) / std::sqrt(2.0 * H6);

        double term1 = exp_safe(u) * erfc_safe(y) / alpha;
        double term2 = exp_safe(v) * erfc_safe(w) / beta;

        return std::sqrt(M_PI / (2.0 * H6)) * (term1 + term2);
      };

      double dt_lo = dt - 0.5 * dt_width;
      double dt_hi = dt + 0.5 * dt_width;

      double f3 = F_tof(dt_hi) - F_tof(dt_lo);

      double result = A_ * (f1 * f2 * f3 / norm_factor);
      if (!std::isfinite(result)) result = 0.0;
      return result;
    }

    double get_norm_factor(Eigen::Matrix3d H, double alpha, double beta) const {
      const std::size_t n = coords.size();
      double sum = 0;
      scitbx::af::versa<double, af::c_grid<3>> out(coords.accessor());
      for (std::size_t c_x = 0; c_x < coords.accessor()[0]; ++c_x) {
        for (std::size_t c_y = 0; c_y < coords.accessor()[1]; ++c_y) {
          for (std::size_t c_z = 0; c_z < coords.accessor()[2]; ++c_z) {
            sum +=
              func(coords(c_x, c_y, c_z), H, alpha, beta, 1.0, 1.0, dt_widths[c_z]);
          }
        }
      }
      return sum;
    }

    double calc_A(const Eigen::Matrix3d& H,
                  double alpha,
                  double beta,
                  double norm_factor) const {
      std::vector<double> P(m);
      int count = 0;
      const double eps = 1e-6;
      for (std::size_t c_x = 0; c_x < coords.accessor()[0]; ++c_x) {
        for (std::size_t c_y = 0; c_y < coords.accessor()[1]; ++c_y) {
          for (std::size_t c_z = 0; c_z < coords.accessor()[2]; ++c_z) {
            P[count] = func(
              coords(c_x, c_y, c_z), H, alpha, beta, 1.0, norm_factor, dt_widths[c_z]);
            count++;
          }
        }
      }

      double num = 0.0, den = 0.0;

      for (size_t i = 0; i < m; ++i) {
        double obs = intensities[i];
        double var = background_variances[i];
        if (!std::isfinite(var) || var <= 0.0) var = std::max(obs, 1e-6);
        double w = 1.0 / var;
        num += w * obs * P[i];
        den += w * P[i] * P[i];
      }

      double A = (den > 0.0) ? (num / den) : 0.0;
      return A;
    }

    int operator()(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) const {
      // Get current params
      Eigen::VectorXd xc = clamp_params(x);
      Eigen::Matrix3d H = build_H_from_L(xc);
      double alpha = std::exp(xc[6]);
      double beta = std::exp(xc[7]);

      Eigen::VectorXd current_params(3);
      current_params << H(0, 0), alpha, beta;

      // Update A and norm factors if params have changed
      if (!cache_valid || (current_params - last_params).cwiseAbs().maxCoeff() > 1e-6) {
        cached_norm = get_norm_factor(H, alpha, beta);
        cached_A = calc_A(H, alpha, beta, cached_norm);
        last_params = current_params;
        cache_valid = true;
      }
      double norm_factor = cached_norm;
      double A = cached_A;

      const std::size_t npts = coords.size();
      fvec.resize(npts);

      std::vector<double> P(m);
      double eps = 1e-8;
      int count = 0;
      for (std::size_t c_x = 0; c_x < coords.accessor()[0]; ++c_x) {
        for (std::size_t c_y = 0; c_y < coords.accessor()[1]; ++c_y) {
          for (std::size_t c_z = 0; c_z < coords.accessor()[2]; ++c_z) {
            double model = func(
              coords(c_x, c_y, c_z), H, alpha, beta, A, norm_factor, dt_widths[c_z]);
            double obs = intensities(c_x, c_y, c_z);
            double var = background_variances(c_x, c_y, c_z);
            if (!std::isfinite(var) || var <= 0.0) {
              var = std::max(obs, 1e-6);
            }
            double sigma = std::sqrt(var);
            double diff = (obs - model) / sigma;

            fvec[count] = std::isfinite(diff) ? diff : 1e6;
            count++;
          }
        }
      }
      return 0;
    }

    int df(const Eigen::VectorXd& x, Eigen::MatrixXd& fjac) const {
      int m = coords.size();
      int n = x.size();
      fjac.resize(m, n);

      Eigen::VectorXd f0(m);
      operator()(x, f0);

      double eps = 1e-5;
      for (int j = 0; j < n; ++j) {
        double delta = eps * std::max(1.0, std::abs(x[j]));
        Eigen::VectorXd xh = x;
        xh[j] += delta;
        Eigen::VectorXd fh(m);
        operator()(xh, fh);
        fjac.col(j) = (fh - f0) / delta;
      }
      return 0;
    }
  };

  class GutmannProfile3D {
  public:
    const scitbx::af::versa<vec3<double>, af::c_grid<3>> coords;
    const scitbx::af::versa<double, af::c_grid<3>> intensities;
    const scitbx::af::versa<double, af::c_grid<3>> background_variances;
    scitbx::af::versa<double, af::c_grid<3>> y_norm;
    int n_restarts;
    double intensity_max;
    boost::optional<GutmannProfileFunctor> functor;

    Eigen::VectorXd params;
    std::array<double, 8> min_bounds;
    std::array<double, 8> max_bounds;

    GutmannProfile3D(
      scitbx::af::const_ref<vec3<double>, af::c_grid<3>> coords_,
      const scitbx::af::versa<double, af::c_grid<3>> intensities_,
      const scitbx::af::versa<double, af::c_grid<3>> background_variances_,
      double alpha,
      double beta,
      const std::array<double, 2> alpha_bounds,
      const std::array<double, 2> beta_bounds,
      int n_restarts_)
        : coords(get_rel_coords(coords_, intensities_)),
          intensities(intensities_),
          background_variances(background_variances_),
          n_restarts(n_restarts_) {
      DIALS_ASSERT(coords.size() == intensities.size());

      const std::size_t n = intensities.size();
      // Ensure no negative values
      intensity_max = 1.0;
      if (!(n == 0)) {
        intensity_max = *std::max_element(intensities.begin(), intensities.end());
        if (intensity_max <= 0.0) intensity_max = 1.0;
      }

      y_norm.resize(intensities.accessor());

      for (std::size_t i = 0; i < n; ++i) {
        double v = intensities[i];
        if (!std::isfinite(v) || v < 0.0) v = 0.0;
        y_norm[i] = v;
      }

      // Compute weighted mean
      Eigen::Vector3d mean = Eigen::Vector3d::Zero();
      double weight_sum = 0.0;

      for (std::size_t c_x = 0; c_x < coords.accessor()[0]; ++c_x) {
        for (std::size_t c_y = 0; c_y < coords.accessor()[1]; ++c_y) {
          for (std::size_t c_z = 0; c_z < coords.accessor()[2]; ++c_z) {
            std::size_t idx =
              (c_x * coords.accessor()[1] + c_y) * coords.accessor()[2] + c_z;
            double w = y_norm[idx];
            vec3<double> coord = coords(c_x, c_y, c_z);
            mean += w * Eigen::Vector3d(coord[0], coord[1], coord[2]);
            weight_sum += w;
          }
        }
      }
      if (weight_sum <= 0.0) weight_sum = 1.0;
      mean /= weight_sum;

      // Compute weighted covariance
      Eigen::Matrix3d C = Eigen::Matrix3d::Zero();
      for (std::size_t c_x = 0; c_x < coords.accessor()[0]; ++c_x) {
        for (std::size_t c_y = 0; c_y < coords.accessor()[1]; ++c_y) {
          for (std::size_t c_z = 0; c_z < coords.accessor()[2]; ++c_z) {
            std::size_t idx =
              (c_x * coords.accessor()[1] + c_y) * coords.accessor()[2] + c_z;
            double w = y_norm[idx];
            vec3<double> coord = coords(c_x, c_y, c_z);
            Eigen::Vector3d v(coord[0], coord[1], coord[2]);
            Eigen::Vector3d d = v - mean;
            C += w * (d * d.transpose());
          }
        }
      }
      C /= weight_sum;
      C += Eigen::Matrix3d::Identity() * 1e-6;

      Eigen::Matrix3d Linv = C.inverse().llt().matrixL();
      double l11 = Linv(0, 0);
      double l21 = Linv(1, 0);
      double l31 = Linv(2, 0);
      double l22 = Linv(1, 1);
      double l32 = Linv(2, 1);
      double l33 = Linv(2, 2);

      params.resize(8);
      params << std::log(l11), l21, l31, std::log(l22), l32, std::log(l33),
        std::log(alpha), std::log(beta);
      // std::cout << "TEST params l11: " << l11 << ", l21: " << l21 << ", l31: " << l31
      //           << ", l22: " << l22 << ", l32: " << l32 << ", l33: " << l33
      //           << ", alpha: " << alpha << " , beta: " << beta << std::endl;

      min_bounds = {std::log(1e-6),
                    -1e2,
                    -1e2,
                    std::log(1e-6),
                    -1e2,
                    std::log(1e-6),
                    std::log(alpha_bounds[0]),
                    std::log(beta_bounds[0])};
      max_bounds = {std::log(1e6),
                    1e2,
                    1e2,
                    std::log(1e6),
                    1e2,
                    std::log(1e6),
                    std::log(alpha_bounds[1]),
                    std::log(beta_bounds[1])};

      // Sanity check
      DIALS_ASSERT(alpha >= alpha_bounds[0] && alpha <= alpha_bounds[1]);
      DIALS_ASSERT(beta >= beta_bounds[0] && beta <= beta_bounds[1]);

      functor.emplace(coords, y_norm, background_variances, min_bounds, max_bounds);
    }

    scitbx::af::versa<vec3<double>, af::c_grid<3>> get_rel_coords(
      scitbx::af::const_ref<vec3<double>, af::c_grid<3>> coords_,
      const scitbx::af::versa<double, af::c_grid<3>> intensities_) {
      double sumI = 0.0;

      double max_intensity = -1.;
      std::size_t coord_idx = 0;
      for (std::size_t i = 0; i < intensities_.size(); ++i) {
        double w = intensities_[i];
        if (!std::isfinite(w) || w <= 0.0) continue;
        if (w > max_intensity) {
          max_intensity = w;
          coord_idx = i;
        }
        sumI += w;
      }
      vec3<double> peak = coords_[coord_idx];
      // Shift coordinates relative to peak
      scitbx::af::versa<vec3<double>, af::c_grid<3>> rel_coords(coords_.accessor());
      for (std::size_t i = 0; i < coords_.size(); ++i) {
        rel_coords[i] = vec3<double>(
          coords_[i][0] - peak[0], coords_[i][1] - peak[1], coords_[i][2] - peak[2]);
      }

      return rel_coords;
    }

    scitbx::af::versa<double, af::c_grid<3>> result() const {
      Eigen::Matrix3d H = functor->build_H_from_L(params);
      double alpha = std::exp(params[6]);
      double beta = std::exp(params[7]);

      double norm_factor = functor->cached_norm;
      double A = functor->cached_A;

      const std::size_t n = coords.size();
      scitbx::af::versa<double, af::c_grid<3>> out(coords.accessor());
      for (std::size_t c_x = 0; c_x < coords.accessor()[0]; ++c_x) {
        for (std::size_t c_y = 0; c_y < coords.accessor()[1]; ++c_y) {
          for (std::size_t c_z = 0; c_z < coords.accessor()[2]; ++c_z) {
            double out_val = functor->func(coords(c_x, c_y, c_z),
                                           H,
                                           alpha,
                                           beta,
                                           A,
                                           norm_factor,
                                           functor->dt_widths[c_z]);
            out(c_x, c_y, c_z) = out_val;
          }
        }
      }
      return out;
    }

    double calc_intensity() const {
      scitbx::af::versa<double, af::c_grid<3>> r = result();

      const auto& acc = r.accessor();
      std::size_t nx = acc[0], ny = acc[1], nz = acc[2];

      double total = 0;
      for (std::size_t ix = 0; ix < nx; ++ix) {
        for (std::size_t iy = 0; iy < ny; ++iy) {
          for (std::size_t iz = 0; iz < nz; ++iz) {
            double v = r(ix, iy, iz);
            total += v;
          }
        }
      }
      return total;
    }

    bool fit(int maxfev = 200, double xtol = 1e-8, double ftol = 1e-8) {
      typedef Eigen::LevenbergMarquardt<GutmannProfileFunctor, double> LM;

      auto run_single_fit = [&](const Eigen::VectorXd& x_init,
                                double& final_error) -> bool {
        Eigen::VectorXd fvec1(functor->intensities.size());
        Eigen::VectorXd x = x_init;
        (*functor)(x, fvec1);
        final_error = fvec1.squaredNorm();
        LM lm(*functor);
        lm.parameters.maxfev = maxfev;
        lm.parameters.xtol = xtol;
        lm.parameters.ftol = ftol;

        int result = lm.minimize(x);
        if (result < 0) return false;

        x = functor->clamp_params(x);
        Eigen::VectorXd fvec(functor->intensities.size());
        (*functor)(x, fvec);
        final_error = fvec.squaredNorm();

        params = x;
        return true;
      };

      // First fit attempt
      Eigen::VectorXd x0 = params;
      double fit_resid = std::numeric_limits<double>::infinity();
      bool success = run_single_fit(x0, fit_resid);
      std::size_t max_profile_index;
      double I_prf, I_var;

      if (success) {
        I_prf = this->calc_intensity();
        if (this->trust_result(fit_resid, I_prf)) {
          return true;
        }
      }

      // Perturb initial params
      std::mt19937 rng(std::random_device{}());
      std::uniform_real_distribution<double> unit_dist(0.0, 1.0);

      for (int i = 0; i < n_restarts; ++i) {
        Eigen::VectorXd x_try(8);

        double scale_factor = 0.5 + unit_dist(rng);
        for (int j = 0; j < 6; ++j) {
          x_try[j] = x0[j] * scale_factor;
        }

        for (int j = 6; j < 8; ++j) {
          double span = max_bounds[j] - min_bounds[j];
          double rand_frac = (unit_dist(rng) - 0.5) * 0.1;
          double perturbed = x0[j] + rand_frac * span;
          x_try[j] = std::max(min_bounds[j], std::min(perturbed, max_bounds[j]));
        }

        success = run_single_fit(x_try, fit_resid);
        if (!success) continue;
        I_prf = this->calc_intensity();
        if (this->trust_result(fit_resid, I_prf)) {
          return true;
        }
      }

      return false;
    }

    bool trust_result(double error, double I_prf) const {
      /*
       * Tests to check the fit is reasonable
       */

      if (!std::isfinite(error) || error <= 0.0) {
        return false;
      }

      // Check reasonable values
      if (I_prf < 1e-7) {
        return false;
      }

      double alpha = std::exp(params[6]);
      double beta = std::exp(params[7]);

      // Check positive H
      Eigen::Matrix3d H = functor->build_H_from_L(params);
      Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eig(H);
      if (eig.eigenvalues().minCoeff() <= 0.0) {
        return false;
      }

      const double exp_limit = 1e50;
      const double mean_limit = 1e20;

      const auto& acc = coords.accessor();
      const std::size_t nx = acc[0];
      const std::size_t ny = acc[1];
      const std::size_t nz = acc[2];

      double profile_peak = 0, data_peak = 0, sum_val = 0;
      double num = 0, denom_y = 0, denom_m = 0;
      double norm_factor = functor->cached_norm;
      double A = functor->cached_A;

      // Check all values sensible
      for (std::size_t ix = 0; ix < nx; ++ix) {
        for (std::size_t iy = 0; iy < ny; ++iy) {
          for (std::size_t iz = 0; iz < nz; ++iz) {
            double data_val = y_norm(ix, iy, iz);
            double val = functor->func(coords(ix, iy, iz),
                                       H,
                                       alpha,
                                       beta,
                                       A,
                                       norm_factor,
                                       functor->dt_widths[iz]);

            if (!std::isfinite(val)) {
              return false;
            }

            data_peak = std::max(data_peak, data_val);
            profile_peak = std::max(profile_peak, val);
            sum_val += val;

            num += data_val * val;
            denom_y += data_val * data_val;
            denom_m += val * val;

            if (val > exp_limit) {
              return false;
            }
          }
        }
      }

      double mean_val = sum_val / coords.size();
      if (mean_val > mean_limit) {
        return false;
      }

      // Check correlation with data
      double corr = num / std::sqrt(denom_y * denom_m + 1e-12);
      if (corr < 0.9) {
        return false;
      }

      // Check peak height is within 10% of data peak
      if (std::abs(profile_peak - data_peak) > data_peak * 0.1) {
        return false;
      }

      return true;
    }
  };

  bool fit_profile3d(
    scitbx::af::const_ref<vec3<double>, af::c_grid<3>> coords,
    const scitbx::af::versa<double, af::c_grid<3>> intensities,
    const scitbx::af::versa<double, af::c_grid<3>> background_variances,
    TOFProfile3DParams& profile_params,
    double& I_prf_out,
    boost::optional<scitbx::af::versa<double, af::c_grid<3>>> profile_3d_out =
      boost::none) {
    // Fit profile
    const std::array<double, 2> alpha_bounds = {profile_params.alpha_min,
                                                profile_params.alpha_max};
    const std::array<double, 2> beta_bounds = {profile_params.beta_min,
                                               profile_params.beta_max};

    GutmannProfile3D profile(coords,
                             intensities,
                             background_variances,
                             profile_params.alpha,
                             profile_params.beta,
                             alpha_bounds,
                             beta_bounds,
                             profile_params.n_restarts);

    bool profile_success = true;
    if (profile_params.optimize_profile) {
      profile_success = profile.fit();
    }

    if (profile_success) {
      double I_prf = profile.calc_intensity();
      I_prf_out = I_prf;
      if (!profile_3d_out) {
        return profile_success;
      }

      scitbx::af::versa<double, af::c_grid<3>> pred = profile.result();
      scitbx::af::versa<double, af::c_grid<3>> profile_3d = *profile_3d_out;

      DIALS_ASSERT(pred.accessor().all_eq(profile_3d.accessor()));
      for (std::size_t x = 0; x < pred.accessor()[0]; ++x) {
        for (std::size_t y = 0; y < pred.accessor()[1]; ++y) {
          for (std::size_t z = 0; z < pred.accessor()[2]; ++z) {
            profile_3d(x, y, z) = pred(x, y, z);
          }
        }
      }
      return profile_success;
    }
    return false;
  }

}}  // namespace dials::algorithms

#endif /* DIALS_ALGORITHMS_INTEGRATION_TOF_GUTMANNPROFILE3D_H */

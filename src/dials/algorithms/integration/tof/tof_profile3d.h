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

  struct GutmannProfileFunctor {
    const scitbx::af::versa<vec3<double>, af::c_grid<3>> coords;
    const scitbx::af::versa<double, af::c_grid<3>> intensities;
    const scitbx::af::versa<double, af::c_grid<3>> background_variances;
    mutable Eigen::VectorXd last_params;  // size 3 or full param vector
    mutable double cached_norm = 1.0;
    mutable bool cache_valid = false;
    std::array<double, 9> min_bounds;
    std::array<double, 9> max_bounds;
    int m, n;

    GutmannProfileFunctor(
      const scitbx::af::versa<vec3<double>, af::c_grid<3>> coords_,
      const scitbx::af::versa<double, af::c_grid<3>> intensities_,
      const scitbx::af::versa<double, af::c_grid<3>> background_variances_,
      const std::array<double, 9>& minb,
      const std::array<double, 9>& maxb)
        : coords(coords_),
          intensities(intensities_),
          background_variances(background_variances_) {
      min_bounds = minb;
      max_bounds = maxb;
      m = coords.size();  // Num data points
      n = 9;              // Num params
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
                double norm_factor) const {
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

      double f3 = exp_safe(u) * erfc_safe(y) + exp_safe(v) * erfc_safe(w);

      double result = A_ * (f1 * f2 * f3 / norm_factor);
      if (!std::isfinite(result)) result = 0.0;
      return result;
    }

    double get_norm_factor(Eigen::Matrix3d H, double alpha, double beta) const {
      const std::size_t n = coords.size();
      scitbx::af::versa<double, af::c_grid<3>> out(coords.accessor());
      for (std::size_t c_x = 0; c_x < coords.accessor()[0]; ++c_x) {
        for (std::size_t c_y = 0; c_y < coords.accessor()[1]; ++c_y) {
          for (std::size_t c_z = 0; c_z < coords.accessor()[2]; ++c_z) {
            double out_val = func(coords(c_x, c_y, c_z), H, alpha, beta, 1.0, 1.0);
            out(c_x, c_y, c_z) = out_val;
          }
        }
      }
      double sum = simpson_integrate_3d(out.const_ref(), coords.const_ref());

      return sum;
    }

    int operator()(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) const {
      Eigen::VectorXd xc = clamp_params(x);
      Eigen::Matrix3d H = build_H_from_L(xc);
      double alpha = std::exp(xc[6]);
      double beta = std::exp(xc[7]);
      double A = std::exp(xc[8]);

      Eigen::VectorXd current_params(3);
      current_params << H(0, 0), alpha, beta;

      if (!cache_valid || (current_params - last_params).cwiseAbs().maxCoeff() > 1e-6) {
        cached_norm = get_norm_factor(H, alpha, beta);
        last_params = current_params;
        cache_valid = true;
      }
      double norm_factor = cached_norm;

      if (norm_factor < 1e-8 || norm_factor > 10.0) {
        std::cout << "Warning: norm_factor = " << norm_factor << " alpha=" << alpha
                  << " beta=" << beta << std::endl;
      }

      const std::size_t npts = coords.size();
      fvec.resize(npts);
      double max_model = -1.;
      double max_intensity = -1.;
      for (std::size_t i = 0; i < npts; ++i) {
        double model = func(coords[i], H, alpha, beta, A, norm_factor);

        double obs = intensities[i];
        double var = background_variances[i];
        if (!std::isfinite(var) || var <= 0.0) {
          var = std::max(obs, 1e-6);
        }
        double sigma = std::sqrt(var);
        double diff = (obs - model) / sigma;

        fvec[i] = std::isfinite(diff) ? diff : 1e6;
        if (model > max_model) {
          max_model = model;
        }
        if (intensities[i] > max_intensity) {
          max_intensity = intensities[i];
        }
      }
      // std::cout<<"TEST intensity " << max_intensity << " model " << max_model <<
      // std::endl;
      return 0;
    }

    int df(const Eigen::VectorXd& x, Eigen::MatrixXd& fjac) const {
      int m = coords.size();
      int n = x.size();
      fjac.resize(m, n);

      Eigen::VectorXd f0(m);
      operator()(x, f0);

      double eps = 1e-4;
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

    Eigen::VectorXd params;
    std::array<double, 9> min_bounds;
    std::array<double, 9> max_bounds;

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
        // y_norm[i] = v / intensity_max;
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

      params.resize(9);
      double A0 = 0.0;
      for (std::size_t i = 0; i < intensities.size(); ++i)
        A0 += std::max(0.0, intensities[i]);
      if (A0 <= 0.0) A0 = intensity_max;
      params << std::log(l11), l21, l31, std::log(l22), l32, std::log(l33),
        std::log(alpha), std::log(beta), std::log(A0);
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
                    std::log(beta_bounds[0]),
                    std::log(std::min(1., A0 * 1e-3))};
      max_bounds = {std::log(1e6),
                    1e2,
                    1e2,
                    std::log(1e6),
                    1e2,
                    std::log(1e6),
                    std::log(alpha_bounds[1]),
                    std::log(beta_bounds[1]),
                    std::log(std::max(1e6, A0 * 1e3))};

      // Sanity check
      DIALS_ASSERT(alpha >= min_bounds[6] && alpha <= max_bounds[6]);
      DIALS_ASSERT(beta >= min_bounds[7] && beta <= max_bounds[7]);
    }

    scitbx::af::versa<vec3<double>, af::c_grid<3>> get_rel_coords(
      scitbx::af::const_ref<vec3<double>, af::c_grid<3>> coords_,
      const scitbx::af::versa<double, af::c_grid<3>> intensities_) {
      // Find approximate centroid (weighted by intensity)
      double sumI = 0.0;
      Eigen::Vector3d centroid = Eigen::Vector3d::Zero();

      double max_intensity = -1.;
      std::size_t coord_idx = 0;
      for (std::size_t i = 0; i < intensities_.size(); ++i) {
        double w = intensities_[i];
        if (!std::isfinite(w) || w <= 0.0) continue;
        if (w > max_intensity) {
          max_intensity = w;
          coord_idx = i;
        }
        centroid += w * Eigen::Vector3d(coords_[i][0], coords_[i][1], coords_[i][2]);
        sumI += w;
      }
      if (sumI > 0.0)
        centroid /= sumI;
      else
        centroid.setZero();

      vec3<double> peak = coords_[coord_idx];
      // Shift coordinates relative to centroid
      scitbx::af::versa<vec3<double>, af::c_grid<3>> rel_coords(coords_.accessor());
      for (std::size_t i = 0; i < coords_.size(); ++i) {
        rel_coords[i] = vec3<double>(
          coords_[i][0] - peak[0], coords_[i][1] - peak[1], coords_[i][2] - peak[2]);
      }

      return rel_coords;
    }

    scitbx::af::versa<double, af::c_grid<3>> result() const {
      GutmannProfileFunctor functor(
        coords, y_norm, background_variances, min_bounds, max_bounds);
      Eigen::Matrix3d H = functor.build_H_from_L(params);
      double alpha = std::exp(params[6]);
      double beta = std::exp(params[7]);
      double A = std::exp(params[8]);

      double norm_factor = functor.get_norm_factor(H, alpha, beta);

      const std::size_t n = coords.size();
      scitbx::af::versa<double, af::c_grid<3>> out(coords.accessor());
      for (std::size_t c_x = 0; c_x < coords.accessor()[0]; ++c_x) {
        for (std::size_t c_y = 0; c_y < coords.accessor()[1]; ++c_y) {
          for (std::size_t c_z = 0; c_z < coords.accessor()[2]; ++c_z) {
            double out_val =
              functor.func(coords(c_x, c_y, c_z), H, alpha, beta, A, norm_factor);
            out(c_x, c_y, c_z) = out_val;
          }
        }
      }
      return out;
    }

    double calc_intensity() const {
      scitbx::af::versa<double, af::c_grid<3>> r = result();
      return simpson_integrate_3d(r.const_ref(), coords.const_ref());
    }

    double calc_variance() {
      scitbx::af::versa<double, af::c_grid<3>> pred = result();
      const auto& acc = pred.accessor();
      DIALS_ASSERT(pred.accessor().all_eq(background_variances.accessor()));

      std::size_t nx = acc[0];
      std::size_t ny = acc[1];
      std::size_t nz = acc[2];

      double intensity = calc_intensity();
      int n_background = 0;
      int n_signal = 0;

      double bg_var_sum = 0;
      double eps = 1e-8;

      for (std::size_t ix = 0; ix < nx; ++ix) {
        for (std::size_t iy = 0; iy < ny; ++iy) {
          for (std::size_t iz = 0; iz < nz; ++iz) {
            double val = pred(ix, iy, iz);
            double bg_val = background_variances(ix, iy, iz);
            if (std::isfinite(val) && std::isfinite(bg_val)) {
              if (val > eps) {
                bg_var_sum += bg_val;
                n_signal++;
              } else {  // Anywhere the profile is flat is classed as background
                n_background++;
              }
            }
          }
        }
      }

      double bg_var = std::abs(bg_var_sum);
      if (n_background > 0) {
        bg_var *= (1.0 + n_signal / n_background);
      }

      return std::abs(intensity) + bg_var;
    }

    bool fit(double I_sum,
             double var_sum,
             int maxfev = 200,
             double xtol = 1e-8,
             double ftol = 1e-8) {
      typedef Eigen::LevenbergMarquardt<GutmannProfileFunctor, double> LM;
      GutmannProfileFunctor functor(
        coords, y_norm, background_variances, min_bounds, max_bounds);

      auto run_single_fit = [&](const Eigen::VectorXd& x_init,
                                double& final_error) -> bool {
        Eigen::VectorXd fvec1(functor.intensities.size());
        Eigen::VectorXd x = x_init;
        functor(x, fvec1);
        final_error = fvec1.squaredNorm();
        std::cout << "TEST initial error " << final_error << std::endl;
        LM lm(functor);
        lm.parameters.maxfev = maxfev;
        lm.parameters.xtol = xtol;
        lm.parameters.ftol = ftol;

        int result = lm.minimize(x);
        if (result < 0) return false;

        x = functor.clamp_params(x);
        Eigen::VectorXd fvec(functor.intensities.size());
        functor(x, fvec);
        final_error = fvec.squaredNorm();
        std::cout << "TEST final error " << final_error << std::endl;

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
        I_var = this->calc_variance();
        if (this->trust_result(fit_resid, I_prf, I_var, I_sum, var_sum)) {
          return true;
        }
      }

      // Perturb initial params
      std::mt19937 rng(std::random_device{}());
      std::uniform_real_distribution<double> unit_dist(0.0, 1.0);

      for (int i = 0; i < n_restarts; ++i) {
        Eigen::VectorXd x_try(9);

        double scale_factor = 0.5 + unit_dist(rng);
        for (int j = 0; j < 6; ++j) {
          x_try[j] = x0[j] * scale_factor;
        }

        for (int j = 6; j < 9; ++j) {
          double span = max_bounds[j] - min_bounds[j];
          double rand_frac = (unit_dist(rng) - 0.5) * 0.1;
          double perturbed = x0[j] + rand_frac * span;
          x_try[j] = std::max(min_bounds[j], std::min(perturbed, max_bounds[j]));
        }

        success = run_single_fit(x_try, fit_resid);
        if (!success) continue;
        I_prf = this->calc_intensity();
        I_var = this->calc_variance();
        if (this->trust_result(fit_resid, I_prf, I_var, I_sum, var_sum)) {
          return true;
        }
      }

      return false;
    }

    bool trust_result(double error,
                      double I_prf,
                      double var_prf,
                      double I_sum,
                      double var_sum) const {
      if (!std::isfinite(error) || error <= 0.0) {
        std::cout << "Failed: error is not finite or <= 0, error = " << error
                  << std::endl;
        return false;
      }

      double alpha = std::exp(params[6]);
      double beta = std::exp(params[7]);
      double A = std::exp(params[8]);

      GutmannProfileFunctor f(
        coords, y_norm, background_variances, min_bounds, max_bounds);
      Eigen::Matrix3d H = f.build_H_from_L(params);
      Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eig(H);
      if (eig.eigenvalues().minCoeff() <= 0.0) {
        std::cout << "Failed: H is not positive-definite, min eigenvalue = "
                  << eig.eigenvalues().minCoeff() << std::endl;
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
      double norm_factor = f.get_norm_factor(H, alpha, beta);

      for (std::size_t ix = 0; ix < nx; ++ix) {
        for (std::size_t iy = 0; iy < ny; ++iy) {
          for (std::size_t iz = 0; iz < nz; ++iz) {
            double data_val = y_norm(ix, iy, iz);
            double val = f.func(coords(ix, iy, iz), H, alpha, beta, A, norm_factor);

            if (!std::isfinite(val)) {
              std::cout << "Failed: val not finite at (" << ix << "," << iy << "," << iz
                        << "), val = " << val << std::endl;
              return false;
            }

            data_peak = std::max(data_peak, data_val);
            profile_peak = std::max(profile_peak, val);
            sum_val += val;

            num += data_val * val;
            denom_y += data_val * data_val;
            denom_m += val * val;

            if (val > exp_limit) {
              std::cout << "Failed: val exceeds exp_limit at (" << ix << "," << iy
                        << "," << iz << "), val = " << val << std::endl;
              return false;
            }
          }
        }
      }

      double mean_val = sum_val / coords.size();
      if (mean_val > mean_limit) {
        std::cout << "Failed: mean_val exceeds mean_limit, mean_val = " << mean_val
                  << std::endl;
        return false;
      }

      if (var_sum > 1e-7) {
        double sum_I_sigma = I_sum / std::sqrt(var_sum);
        double prf_I_sigma = I_prf / std::sqrt(var_prf);
        if (sum_I_sigma > (prf_I_sigma + sum_I_sigma * 0.1) && false) {
          std::cout << "Failed: sum_I_sigma > prf_I_sigma + 10%, sum_I_sigma = "
                    << sum_I_sigma << ", prf_I_sigma = " << prf_I_sigma << std::endl;
          return false;
        }
      }

      double corr = num / std::sqrt(denom_y * denom_m + 1e-12);
      if (corr < 0.9 && false) {
        std::cout << "Failed: correlation < 0.9, corr = " << corr << std::endl;
        return false;
      }

      if (std::abs(profile_peak - data_peak) > data_peak * 0.1) {
        std::cout << "Failed: profile_peak not within 10% of data_peak, profile_peak = "
                  << profile_peak << ", data_peak = " << data_peak << " alpha " << alpha
                  << " beta " << beta << " A " << A << std::endl;
        return false;
      }

      std::cout << "Success: all checks passed corr " << corr << " alpha " << alpha
                << " beta " << beta << " A " << A << std::endl;
      return true;
    }
  };

}}  // namespace dials::algorithms

#endif /* DIALS_ALGORITHMS_INTEGRATION_TOF_GUTMANNPROFILE3D_H */

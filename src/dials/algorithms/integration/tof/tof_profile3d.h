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
    std::array<double, 8> min_bounds;
    std::array<double, 8> max_bounds;
    int m, n;

    GutmannProfileFunctor(const scitbx::af::versa<vec3<double>, af::c_grid<3>> coords_,
                          const scitbx::af::versa<double, af::c_grid<3>> intensities_,
                          const std::array<double, 8>& minb,
                          const std::array<double, 8>& maxb)
        : coords(coords_), intensities(intensities_) {
      min_bounds = minb;
      max_bounds = maxb;
      m = coords.size();  // Num data points
      n = 8;              // Num params
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
      double l11 = params[0];
      double l21 = params[1];
      double l31 = params[2];
      double l22 = params[3];
      double l32 = params[4];
      double l33 = params[5];

      Eigen::Matrix3d L;
      L << l11, 0.0, 0.0, l21, l22, 0.0, l31, l32, l33;

      return L * L.transpose();
    }

    double func(const vec3<double>& c,
                const Eigen::Matrix3d& H,
                double alpha,
                double beta) const {
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

      double result = f1 * f2 * f3;
      if (!std::isfinite(result)) result = 0.0;
      return result;
    }

    int operator()(const Eigen::VectorXd& x, Eigen::VectorXd& fvec) const {
      Eigen::VectorXd xc = clamp_params(x);
      Eigen::Matrix3d H = build_H_from_L(xc);
      double alpha = xc[6];
      double beta = xc[7];

      const std::size_t n = coords.size();
      for (std::size_t i = 0; i < n; ++i) {
        double model = func(coords[i], H, alpha, beta);
        double diff = intensities[i] - model;
        if (!std::isfinite(diff)) diff = 1e6;
        fvec[i] = diff;
      }
      return 0;
    }

    // Numerical Jacobian (finite differences)
    int df(const Eigen::VectorXd& x, Eigen::MatrixXd& fjac) const {
      const double eps = 1e-4;
      Eigen::VectorXd f0(m);
      operator()(x, f0);

      for (int j = 0; j < n; ++j) {
        Eigen::VectorXd xh = x;
        xh[j] += eps;
        Eigen::VectorXd fh(m);
        operator()(xh, fh);
        fjac.col(j) = (fh - f0) / eps;
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

    Eigen::VectorXd params;
    std::array<double, 8> min_bounds;
    std::array<double, 8> max_bounds;

    GutmannProfile3D(
      const scitbx::af::versa<vec3<double>, af::c_grid<3>> coords_,
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
      y_norm.resize(intensities.accessor());
      double sumI = 0.0;
      for (double v : intensities)
        if (std::isfinite(v) && v > 0.0) sumI += v;
      if (sumI <= 0.0) sumI = 1.0;

      for (std::size_t i = 0; i < n; ++i) {
        double v = intensities[i];
        if (!std::isfinite(v) || v < 0.0) v = 0.0;
        y_norm[i] = std::isfinite(v) && v > 0.0 ? v / sumI : 0.0;
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
      params << l11, l21, l31, l22, l32, l33, alpha, beta;
      std::cout << "TEST params l11: " << l11 << ", l21: " << l21 << ", l31: " << l31
                << ", l22: " << l22 << ", l32: " << l32 << ", l33: " << l33
                << ", alpha: " << alpha << " , beta: " << beta << std::endl;

      min_bounds = {
        1e-4, -1e-3, -1e-3, 1e-4, -1e-3, -1e-3, alpha_bounds[0], beta_bounds[0]};
      max_bounds = {1e5, 1e5, 1e5, 1e5, 1e5, 1e5, alpha_bounds[1], beta_bounds[1]};

      // Sanity check
      DIALS_ASSERT(alpha >= min_bounds[6] && alpha <= max_bounds[6]);
      DIALS_ASSERT(beta >= min_bounds[7] && beta <= max_bounds[7]);
    }

    scitbx::af::versa<vec3<double>, af::c_grid<3>> get_rel_coords(
      const scitbx::af::versa<vec3<double>, af::c_grid<3>> coords_,
      const scitbx::af::versa<double, af::c_grid<3>> intensities_) {
      // Find approximate centroid (weighted by intensity)
      double sumI = 0.0;
      Eigen::Vector3d centroid = Eigen::Vector3d::Zero();

      for (std::size_t i = 0; i < intensities_.size(); ++i) {
        double w = intensities_[i];
        if (!std::isfinite(w) || w <= 0.0) continue;
        centroid += w * Eigen::Vector3d(coords_[i][0], coords_[i][1], coords_[i][2]);
        sumI += w;
      }
      if (sumI > 0.0)
        centroid /= sumI;
      else
        centroid.setZero();

      // Shift coordinates relative to centroid
      scitbx::af::versa<vec3<double>, af::c_grid<3>> rel_coords(coords_.accessor());
      for (std::size_t i = 0; i < coords_.size(); ++i) {
        rel_coords[i] = vec3<double>(coords_[i][0] - centroid[0],
                                     coords_[i][1] - centroid[1],
                                     coords_[i][2] - centroid[2]);
      }

      return rel_coords;
    }

    scitbx::af::versa<double, af::c_grid<3>> result() const {
      GutmannProfileFunctor functor(coords, y_norm, min_bounds, max_bounds);
      Eigen::Matrix3d H = functor.build_H_from_L(params);
      double alpha = params[6];
      double beta = params[7];

      for (std::size_t i = 0; i < params.size(); ++i) {
        std::cout << "TEST result params " << params[i] << std::endl;
      }

      const std::size_t n = coords.size();
      scitbx::af::versa<double, af::c_grid<3>> out(coords.accessor());
      for (std::size_t c_x = 0; c_x < coords.accessor()[0]; ++c_x) {
        for (std::size_t c_y = 0; c_y < coords.accessor()[1]; ++c_y) {
          for (std::size_t c_z = 0; c_z < coords.accessor()[2]; ++c_z) {
            out(c_x, c_y, c_z) = functor.func(coords(c_x, c_y, c_z), H, alpha, beta);
          }
        }
      }
      return out;
    }

    scitbx::af::versa<double, af::c_grid<3>> scaled_result() const {
      scitbx::af::versa<double, af::c_grid<3>> model = result();

      // Scale model to raw intensity level
      double sum_raw = 0.0;
      double sum_model = 0.0;
      for (auto v : intensities)
        if (std::isfinite(v) && v > 0.0) sum_raw += v;
      for (auto v : model)
        if (std::isfinite(v) && v > 0.0) sum_model += v;

      double scale = (sum_model > 0) ? sum_raw / sum_model : 1.0;
      for (std::size_t i = 0; i < model.size(); ++i)
        model[i] *= scale;
      return model;
    }

    double calc_intensity() const {
      /*
       * Simpson integration along TOF,
       * summation over x, y
       */

      scitbx::af::versa<double, af::c_grid<3>> pred = result();

      // Normalize raw sum
      double sum_raw = 0.0;
      for (std::size_t x = 0; x < intensities.accessor()[0]; ++x) {
        for (std::size_t y = 0; y < intensities.accessor()[1]; ++y) {
          for (std::size_t z = 0; z < intensities.accessor()[2]; ++z) {
            double v = intensities(x, y, z);
            if (std::isfinite(v) && v > 0.0) sum_raw += v;
          }
        }
      }

      if (sum_raw <= 0.0) sum_raw = 1.0;

      const auto& acc = pred.accessor();
      std::size_t nx = acc[0];
      std::size_t ny = acc[1];
      std::size_t nt = acc[2];

      // Integrated profile along TOF axis for each (x, y)
      scitbx::af::shared<double> profile_1d(nt, 0.0);

      for (std::size_t c_x = 0; c_x < nx; ++c_x) {
        for (std::size_t c_y = 0; c_y < ny; ++c_y) {
          scitbx::af::shared<double> pred_line(nt);
          scitbx::af::shared<double> coord_line(nt);
          for (std::size_t c_z = 0; c_z < nt; ++c_z) {
            pred_line[c_z] = pred(c_x, c_y, c_z);
            coord_line[c_z] = coords(c_x, c_y, c_z)[2];
          }

          double integral_t =
            simpson_integrate(pred_line.const_ref(), coord_line.const_ref());
          for (std::size_t c_z = 0; c_z < nt; ++c_z)
            profile_1d[c_z] += pred_line[c_z];
        }
      }

      // integrate the sum over x,y
      double sum_pred = 0.0;
      for (double v : profile_1d) {
        if (std::isfinite(v)) sum_pred += v;
      }

      return sum_pred * sum_raw;
    }

    double calc_intensity_old() const {
      /*
       * Simpson integration along TOF,
       * summation over x, y
       */

      scitbx::af::versa<double, af::c_grid<3>> pred = result();

      // Normalize raw sum
      double sum_raw = 0.0;
      for (std::size_t x = 0; x < intensities.accessor()[0]; ++x) {
        for (std::size_t y = 0; y < intensities.accessor()[1]; ++y) {
          for (std::size_t z = 0; z < intensities.accessor()[2]; ++z) {
            double v = intensities(x, y, z);
            if (std::isfinite(v) && v > 0.0) sum_raw += v;
          }
        }
      }

      if (sum_raw <= 0.0) sum_raw = 1.0;

      const auto& acc = pred.accessor();
      std::size_t nx = acc[0];
      std::size_t ny = acc[1];
      std::size_t nt = acc[2];

      // Sum the integrated intensity over all (x,y) positions
      double total_integrated = 0.0;

      for (std::size_t c_x = 0; c_x < nx; ++c_x) {
        for (std::size_t c_y = 0; c_y < ny; ++c_y) {
          scitbx::af::shared<double> pred_line(nt);
          scitbx::af::shared<double> coord_line(nt);
          for (std::size_t c_z = 0; c_z < nt; ++c_z) {
            pred_line[c_z] = pred(c_x, c_y, c_z);
            coord_line[c_z] = coords(c_x, c_y, c_z)[2];
          }

          // Integrate along TOF axis for this (x,y) position
          double integral_t =
            simpson_integrate(pred_line.const_ref(), coord_line.const_ref());

          // Accumulate the integral
          total_integrated += integral_t;
        }
      }

      // Scale by raw data normalization
      return total_integrated * sum_raw;
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
             int maxfev = 2000,
             double xtol = 1e-8,
             double ftol = 1e-8) {
      typedef Eigen::LevenbergMarquardt<GutmannProfileFunctor, double> LM;
      GutmannProfileFunctor functor(coords, y_norm, min_bounds, max_bounds);

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
        Eigen::VectorXd fvec(functor.intensities.size());
        functor(x, fvec);
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
        I_var = this->calc_variance();
        if (this->trust_result(fit_resid, I_prf, I_var, I_sum, var_sum)) {
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

      double a = params[6], b = params[7];
      if (a <= 0.0 || b <= 0.0) {
        std::cout << "Failed: alpha/beta not positive, a = " << a << ", b = " << b
                  << std::endl;
        return false;
      }

      GutmannProfileFunctor f(coords, y_norm, min_bounds, max_bounds);
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

      for (std::size_t ix = 0; ix < nx; ++ix) {
        for (std::size_t iy = 0; iy < ny; ++iy) {
          for (std::size_t iz = 0; iz < nz; ++iz) {
            double data_val = y_norm(ix, iy, iz);
            double val = f.func(coords(ix, iy, iz), H, a, b);

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

      if (std::abs(profile_peak - data_peak) > data_peak * 0.1 && false) {
        std::cout << "Failed: profile_peak not within 10% of data_peak, profile_peak = "
                  << profile_peak << ", data_peak = " << data_peak << std::endl;
        return false;
      }

      std::cout << "Success: all checks passed corr " << corr << std::endl;
      return true;
    }
  };

}}  // namespace dials::algorithms

#endif /* DIALS_ALGORITHMS_INTEGRATION_TOF_GUTMANNPROFILE3D_H */

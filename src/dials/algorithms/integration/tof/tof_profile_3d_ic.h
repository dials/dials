#ifndef DIALS_ALGORITHMS_INTEGRATION_TOF_PROFILE_3D_IC_H
#define DIALS_ALGORITHMS_INTEGRATION_TOF_PROFILE_3D_IC_H

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
Ikeda-Carpenter (IC) (ToF) × Bivariate Gaussian (BVG) model (XY)

Includes a discrete convolution  to map the raw ic profile to match instrument bin
resolution

based on
Ikeda, S.I Carpenter J.M, Wide-energy-range, high-resolution measurements of
neutron pulse shapes of polyethylene moderators, Nucl. Instrum. Methods Phys. Res. A.
239(3), 536-544, (1985) https://doi.org/10.1016/0168-9002(85)90033-6


Mantid's IkedaCarpenterPV method:
https://docs.mantidproject.org/nightly/fitting/fitfunctions/IkedaCarpenterPV.html


*/

namespace dials { namespace algorithms {

  /*
   * Holds params required for TOFProfile3DIC
   */

  struct TOFProfile3DICParams {
    // Ikeda-Carpenter parameters
    double A, A_min, A_max;  // fast-neutron decay constant
    double B, B_min, B_max;  // slow-neutron decay constant
    double R, R_min, R_max;  // fast/slow neutron ratio
                             //
    // BVG parameters
    double SigX_min, SigX_max;
    double SigY_min, SigY_max;
    double SigP, SigP_min, SigP_max;  // spatial correlation
                                      //
    // Convolution parameters (optimised by default but can be fixed with
    // optimize_convolution_params=false)
    double HatWidth;  // half-width of rectangular convolution kernel in ToF bins
    double KConv;     // Gaussian kernel decay rate per ToF bin^2 (sigma =
                      // 1/sqrt(2*KConv) bins)

    // Control
    int n_restarts;         // number of attempts when fitting
    bool optimize_profile;  // If false the profile is generated with input params
    bool optimize_convolution_params;  // If false, fix HatWidth and KConv
    bool optimize_moderator_params;    // If false, fix the Ikeda-Carpenter A, B, R
                                       // (moderator) params so only the per-reflection
                                       // shape/position params are fitted
    double max_drift_factor;           // Scales the per-reflection dXdt/dYdt bounds
                                       // relative to the shoebox's spatial extent
                                       // divided by its ToF extent
    bool show_profile_failures;        // Prints debugging information

    TOFProfile3DICParams(double A_,
                         double A_min_,
                         double A_max_,
                         double B_,
                         double B_min_,
                         double B_max_,
                         double R_,
                         double R_min_,
                         double R_max_,
                         double SigX_min_,
                         double SigX_max_,
                         double SigY_min_,
                         double SigY_max_,
                         double SigP_,
                         double SigP_min_,
                         double SigP_max_,
                         double HatWidth_,
                         double KConv_,
                         int n_restarts_,
                         bool optimize_profile_,
                         bool optimize_convolution_params_,
                         bool optimize_moderator_params_,
                         double max_drift_factor_,
                         bool show_profile_failures_)
        : A(A_),
          A_min(A_min_),
          A_max(A_max_),
          B(B_),
          B_min(B_min_),
          B_max(B_max_),
          R(R_),
          R_min(R_min_),
          R_max(R_max_),
          SigX_min(SigX_min_),
          SigX_max(SigX_max_),
          SigY_min(SigY_min_),
          SigY_max(SigY_max_),
          SigP(SigP_),
          SigP_min(SigP_min_),
          SigP_max(SigP_max_),
          HatWidth(HatWidth_),
          KConv(KConv_),
          n_restarts(n_restarts_),
          optimize_profile(optimize_profile_),
          optimize_convolution_params(optimize_convolution_params_),
          optimize_moderator_params(optimize_moderator_params_),
          max_drift_factor(max_drift_factor_),
          show_profile_failures(show_profile_failures_) {}
  };

  template <typename T>
  static scitbx::af::shared<T> ic_raw_t(scitbx::af::const_ref<double> tof,
                                        T A,
                                        T B,
                                        T R,
                                        T T0) {
    /*
     * Ikeda-Carpenter function with unit amplitude for each tof[i].
     */

    const int n = tof.size();
    scitbx::af::shared<T> out(n, T(0.0));
    T dAB = A - B;
    if (std::abs(as_real(dAB)) < 1e-8) dAB = T(as_real(dAB) >= 0.0 ? 1e-8 : -1e-8);
    const T dAB3 = dAB * dAB * dAB;
    const T coeff2 = T(2.0) * R * A * A * B / dAB3;
    for (int i = 0; i < n; ++i) {
      const T dt = T(tof[i]) - T0;
      if (as_real(dt) <= 0.0) continue;
      const T Adt = A * dt;
      const T Bdt = B * dt;
      const T dABdt = dAB * dt;
      const T expA = exp_safe_t(-Adt);
      const T expB = exp_safe_t(-Bdt);
      const T term1 = (T(1.0) - R) * Adt * Adt * expA;
      const T term2 =
        coeff2 * (expB - expA * (T(1.0) + dABdt + T(0.5) * dABdt * dABdt));
      const T val = term1 + term2;
      out[i] = (std::isfinite(as_real(val)) && as_real(val) > 0.0) ? val : T(0.0);
    }
    return out;
  }

  static scitbx::af::shared<double> ic_raw(scitbx::af::const_ref<double> tof,
                                           double A,
                                           double B,
                                           double R,
                                           double T0) {
    return ic_raw_t<double>(tof, A, B, R, T0);
  }

  static scitbx::af::shared<double> make_hat_kernel(int n, double hw) {
    // Normalised hat kernel of length n, half-width hw bins

    scitbx::af::shared<double> k(n, 0.0);
    const int mid = n / 2;
    const int lo = std::max(0, static_cast<int>(std::floor(mid - std::abs(hw))));
    const int hi = std::min(n, static_cast<int>(std::ceil(mid + std::abs(hw))));
    for (int i = lo; i < hi; ++i)
      k[i] = 1.0;
    double sum = 0.0;
    for (int i = 0; i < n; ++i)
      sum += k[i];
    if (sum > 0.0)
      for (int i = 0; i < n; ++i)
        k[i] /= sum;
    else
      k[mid] = 1.0;  // identity kernel when hw=0
    return k;
  }

  static scitbx::af::shared<double> make_gauss_kernel(int n, double kconv) {
    /*
     * Normalised Gaussian kernel of length n.
     * kconv is a decay rate inToF-bin units
     */

    scitbx::af::shared<double> k(n, 0.0);
    const double mid = (n > 1) ? 0.5 * (n - 1) : 0.0;
    double sum = 0.0;
    for (int i = 0; i < n; ++i) {
      const double d = i - mid;
      k[i] = std::exp(-kconv * d * d);
      sum += k[i];
    }
    if (sum > 0.0)
      for (int i = 0; i < n; ++i)
        k[i] /= sum;
    return k;
  }

  static scitbx::af::shared<double> convolve(scitbx::af::const_ref<double> f,
                                             scitbx::af::const_ref<double> kernel) {
    // Discrete convolution

    const int n = f.size();
    const int m = kernel.size();
    const int first = n / 2;
    scitbx::af::shared<double> out(n, 0.0);
    for (int i = 0; i < n; ++i) {
      for (int k = 0; k < n; ++k) {
        const int ki = i + first - k;
        if (ki >= 0 && ki < m) out[i] += f[k] * kernel[ki];
      }
    }
    return out;
  }

  static scitbx::af::shared<double> ic_profile(scitbx::af::const_ref<double> tof,
                                               double A,
                                               double B,
                                               double R,
                                               double T0,
                                               double HatWidth,
                                               double KConv) {
    const int n = tof.size();
    scitbx::af::shared<double> f = ic_raw(tof, A, B, R, T0);
    scitbx::af::shared<double> hat = make_hat_kernel(n, HatWidth);
    scitbx::af::shared<double> gauss = make_gauss_kernel(n, KConv);
    scitbx::af::shared<double> f_hat = convolve(f.const_ref(), hat.const_ref());
    return convolve(f_hat.const_ref(), gauss.const_ref());
  }

  struct Profile3DICFunctor {
    const scitbx::af::versa<vec3<double>, af::c_grid<3>>
      coords;  // (x,y,t) centred at peak
    const scitbx::af::versa<double, af::c_grid<3>> intensities;
    const scitbx::af::versa<double, af::c_grid<3>> background_variances;
    scitbx::af::shared<double> tof_axis;  // 1D ToF, length = accessor()[2]
    double fixed_hat_width, fixed_kconv;
    std::array<double, 11> min_bounds, max_bounds;
    int num_data_points, num_params;
    bool optimize_convolution_params;
    bool optimize_moderator_params;
    mutable double cached_Scale;
    mutable Eigen::VectorXd last_params;
    mutable bool cache_valid;

    /*
     * Parameter mapping:
     *   x[0] = log(A)       fast-neutron decay
     *   x[1] = log(B)       slow-neutron decay
     *   x[2] = R            fast/slow ratio [0,1], clamped
     *   x[3] = T0           time onset, relative to centred peak
     *   x[4] = log(SigX)    BVG sigma x
     *   x[5] = log(SigY)    BVG sigma y
     *   x[6] = SigP         BVG correlation (-1,1), clamped
     *   x[7] = dXdt         BVG centre drift in x per unit ToF (peak walk)
     *   x[8] = dYdt         BVG centre drift in y per unit ToF (peak walk)
     * When optimize_convolution_params:
     *   x[9]  = log(HatWidth)
     *   x[10] = log(KConv)
     */

    Profile3DICFunctor(
      const scitbx::af::versa<vec3<double>, af::c_grid<3>> coords_,
      const scitbx::af::versa<double, af::c_grid<3>> intensities_,
      const scitbx::af::versa<double, af::c_grid<3>> background_variances_,
      const std::array<double, 11>& minb,
      const std::array<double, 11>& maxb,
      double hat_width,
      double kconv,
      bool opt_conv_params,
      bool opt_moderator_params)
        : coords(coords_),
          intensities(intensities_),
          background_variances(background_variances_),
          fixed_hat_width(hat_width),
          fixed_kconv(kconv),
          optimize_convolution_params(opt_conv_params),
          optimize_moderator_params(opt_moderator_params),
          cached_Scale(1.0),
          cache_valid(false) {
      min_bounds = minb;
      max_bounds = maxb;
      num_data_points = coords.size();
      num_params = opt_conv_params ? 11 : 9;
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
                        double& A,
                        double& B,
                        double& R,
                        double& T0,
                        double& SigX,
                        double& SigY,
                        double& SigP,
                        double& dXdt,
                        double& dYdt,
                        double& HatWidth,
                        double& KConv) const {
      // Extracts params from xc

      A = std::exp(xc[0]);
      B = std::exp(xc[1]);
      R = xc[2];
      T0 = xc[3];
      SigX = std::exp(xc[4]);
      SigY = std::exp(xc[5]);
      SigP = xc[6];
      dXdt = xc[7];
      dYdt = xc[8];
      HatWidth = optimize_convolution_params ? std::exp(xc[9]) : fixed_hat_width;
      KConv = optimize_convolution_params ? std::exp(xc[10]) : fixed_kconv;
    }

    double calc_scale(const scitbx::af::shared<double>& ic_conv,
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
            const double p = ic_conv[iz] * bvg;
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
      double A, B, R, T0, SigX, SigY, SigP, dXdt, dYdt, HatWidth, KConv;
      extract_params(xc, A, B, R, T0, SigX, SigY, SigP, dXdt, dYdt, HatWidth, KConv);

      const scitbx::af::shared<double> ic_conv =
        ic_profile(tof_axis.const_ref(), A, B, R, T0, HatWidth, KConv);

      const int nx = coords.accessor()[0];
      const int ny = coords.accessor()[1];
      const int nz = coords.accessor()[2];

      if (!cache_valid || (xc - last_params).cwiseAbs().maxCoeff() > 1e-6) {
        cached_Scale = calc_scale(ic_conv, SigX, SigY, SigP, dXdt, dYdt);
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
            const double model = Scale * ic_conv[iz] * bvg;
            const double obs = intensities(ix, iy, iz);
            double var_b = background_variances(ix, iy, iz);
            if (!std::isfinite(var_b) || var_b < 0.0) var_b = 0.0;
            const double var = std::max(std::abs(obs) + var_b, 1.0);
            const double sigma = std::sqrt(var);
            const double diff = (obs - model) / sigma;
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
        // Hold the Ikeda-Carpenter moderator params (A=0, B=1, R=2) fixed by
        // giving them a zero Jacobian column, so LM leaves them unchanged
        if (!optimize_moderator_params && j <= 2) {
          J.col(j).setZero();
          continue;
        }
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

  class TOFProfile3DIC {
  public:
    const scitbx::af::versa<vec3<double>, af::c_grid<3>> coords;
    const scitbx::af::versa<double, af::c_grid<3>> intensities;
    const scitbx::af::versa<double, af::c_grid<3>> background_variances;
    scitbx::af::versa<double, af::c_grid<3>> y_norm;
    int n_restarts;
    Eigen::VectorXd params;
    std::array<double, 11> min_bounds, max_bounds;
    double tof_range;
    boost::optional<Profile3DICFunctor> functor;

    TOFProfile3DIC(scitbx::af::const_ref<vec3<double>, af::c_grid<3>> coords_,
                   const scitbx::af::versa<double, af::c_grid<3>> intensities_,
                   const scitbx::af::versa<double, af::c_grid<3>> background_variances_,
                   double A,
                   double B,
                   double R,
                   double SigP,
                   double HatWidth,
                   double KConv,
                   const std::array<double, 2>& A_bounds,
                   const std::array<double, 2>& B_bounds,
                   const std::array<double, 2>& R_bounds,
                   const std::array<double, 2>& SigX_bounds,
                   const std::array<double, 2>& SigY_bounds,
                   const std::array<double, 2>& SigP_bounds,
                   int n_restarts_,
                   bool optimize_convolution_params,
                   bool optimize_moderator_params,
                   double max_drift_factor)
        : coords(get_rel_coords(coords_, intensities_)),
          intensities(intensities_),
          background_variances(background_variances_),
          n_restarts(n_restarts_) {
      const int n = intensities.size();
      y_norm.resize(intensities.accessor());
      for (int i = 0; i < n; ++i) {
        const double v = intensities[i];
        y_norm[i] = (std::isfinite(v) && v > 0.0) ? v : 0.0;
      }

      // Estimate BVG sigmas from intensity-weighted spatial variance
      double wsum = 0.0, wx2 = 0.0, wy2 = 0.0;
      for (std::size_t ix = 0; ix < coords.accessor()[0]; ++ix)
        for (std::size_t iy = 0; iy < coords.accessor()[1]; ++iy)
          for (std::size_t iz = 0; iz < coords.accessor()[2]; ++iz) {
            const double w = y_norm(ix, iy, iz);
            const double dx = coords(ix, iy, iz)[0];
            const double dy = coords(ix, iy, iz)[1];
            wsum += w;
            wx2 += w * dx * dx;
            wy2 += w * dy * dy;
          }
      if (wsum <= 0.0) wsum = 1.0;
      const double SigX_est =
        std::sqrt(std::max(wx2 / wsum, SigX_bounds[0] * SigX_bounds[0]));
      const double SigY_est =
        std::sqrt(std::max(wy2 / wsum, SigY_bounds[0] * SigY_bounds[0]));
      const double SigX_init =
        std::min(std::max(SigX_est, SigX_bounds[0]), SigX_bounds[1]);
      const double SigY_init =
        std::min(std::max(SigY_est, SigY_bounds[0]), SigY_bounds[1]);
      const double SigP_init = std::min(std::max(SigP, SigP_bounds[0]), SigP_bounds[1]);

      // T0 init: IC function peaks at T0 + 2/A so centering places the data peak at t=0
      // Set T0 = -2/A so the IC peak aligns with t=0.  T0 must be negative (before
      // the peak), but need not be before tof_min
      // The IC function handles t<=T0 by returning 0, so T0 can lie inside the data
      // window
      const int nz = coords.accessor()[2];
      const double tof_min = coords(0, 0, 0)[2];
      const double tof_max = coords(0, 0, nz - 1)[2];
      tof_range = std::abs(tof_max - tof_min);
      const double T0_min = tof_min - 2.0 * tof_range;
      const double T0_max = -1e-6;  // T0 must be before the peak at t=0
      const double T0_init = std::min(std::max(-2.0 / A, T0_min), T0_max);

      // Drift bounds: cap how far the BVG centre can move over the shoebox's
      // ToF range, relative to the shoebox's own spatial extent (falling back
      // to the SigX/SigY upper bound if the shoebox is degenerate in x or y),
      // scaled by max_drift_factor. dXdt/dYdt default to 0 (no drift).
      const std::size_t nx = coords.accessor()[0];
      const std::size_t ny = coords.accessor()[1];
      const double x_extent = std::abs(coords(nx - 1, 0, 0)[0] - coords(0, 0, 0)[0]);
      const double y_extent = std::abs(coords(0, ny - 1, 0)[1] - coords(0, 0, 0)[1]);
      const double drift_denom = std::max(tof_range, 1e-6);
      const double dXdt_max =
        max_drift_factor * std::max(x_extent, SigX_bounds[1]) / drift_denom;
      const double dYdt_max =
        max_drift_factor * std::max(y_extent, SigY_bounds[1]) / drift_denom;

      // Bounds array: (logA, logB, R, T0, logSigX, logSigY, SigP, dXdt, dYdt,
      // logHatW, logKConv)
      min_bounds = {std::log(A_bounds[0]),
                    std::log(B_bounds[0]),
                    R_bounds[0],
                    T0_min,
                    std::log(SigX_bounds[0]),
                    std::log(SigY_bounds[0]),
                    SigP_bounds[0],
                    -dXdt_max,
                    -dYdt_max,
                    std::log(0.5),
                    std::log(0.02)};
      max_bounds = {std::log(A_bounds[1]),
                    std::log(B_bounds[1]),
                    R_bounds[1],
                    T0_max,
                    std::log(SigX_bounds[1]),
                    std::log(SigY_bounds[1]),
                    SigP_bounds[1],
                    dXdt_max,
                    dYdt_max,
                    std::log(20.0),
                    std::log(5.0)};

      params.resize(optimize_convolution_params ? 11 : 9);
      params[0] = std::log(A);
      params[1] = std::log(B);
      params[2] = std::min(std::max(R, R_bounds[0]), R_bounds[1]);
      params[3] = T0_init;
      params[4] = std::log(SigX_init);
      params[5] = std::log(SigY_init);
      params[6] = SigP_init;
      params[7] = 0.0;  // dXdt: no drift assumed unless the fit finds evidence
      params[8] = 0.0;  // dYdt
      if (optimize_convolution_params) {
        params[9] = std::log(std::max(HatWidth, 0.5));
        params[10] = std::log(std::min(std::max(KConv, 0.02), 5.0));
      }

      functor.emplace(coords,
                      y_norm,
                      background_variances,
                      min_bounds,
                      max_bounds,
                      HatWidth,
                      KConv,
                      optimize_convolution_params,
                      optimize_moderator_params);
    }

    scitbx::af::versa<vec3<double>, af::c_grid<3>> get_rel_coords(
      scitbx::af::const_ref<vec3<double>, af::c_grid<3>> coords_,
      const scitbx::af::versa<double, af::c_grid<3>> intensities_) {
      double max_intensity = -1.0;
      std::size_t peak_idx = 0;
      for (std::size_t i = 0; i < intensities_.size(); ++i) {
        const double w = intensities_[i];
        if (std::isfinite(w) && w > max_intensity) {
          max_intensity = w;
          peak_idx = i;
        }
      }
      const vec3<double> peak = coords_[peak_idx];
      scitbx::af::versa<vec3<double>, af::c_grid<3>> rel(coords_.accessor());
      for (std::size_t i = 0; i < coords_.size(); ++i)
        rel[i] = vec3<double>(
          coords_[i][0] - peak[0], coords_[i][1] - peak[1], coords_[i][2] - peak[2]);
      return rel;
    }

    double get_A() const {
      return std::exp(params[0]);
    }
    double get_B() const {
      return std::exp(params[1]);
    }
    double get_R() const {
      return params[2];
    }
    double get_T0() const {
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
      const double A = get_A(), B = get_B(), R = get_R(), T0 = get_T0();
      const double SigX = get_SigX(), SigY = get_SigY(), SigP = get_SigP();
      const double dXdt = get_dXdt(), dYdt = get_dYdt();
      double HatWidth = functor->fixed_hat_width;
      double KConv = functor->fixed_kconv;
      if (functor->optimize_convolution_params && params.size() >= 11) {
        HatWidth = std::exp(params[9]);
        KConv = std::exp(params[10]);
      }
      const scitbx::af::shared<double> ic_conv =
        ic_profile(functor->tof_axis.const_ref(), A, B, R, T0, HatWidth, KConv);
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
            const double val = Scale * ic_conv[iz] * bvg;
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
      typedef Eigen::LevenbergMarquardt<Profile3DICFunctor, double> LM;

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

      // Log-scale params
      const bool opt_mod = functor->optimize_moderator_params;
      const int n_log = functor->optimize_convolution_params ? 6 : 4;
      const int log_idx[6] = {0, 1, 4, 5, 9, 10};
      auto is_fixed = [&](int j) {
        if (!opt_mod && j <= 2) return true;  // A, B, R fixed
        return false;
      };

      // Params that are not log scale
      auto perturb_position = [&](Eigen::VectorXd& x_try, double scale) {
        for (int j : {3, 7, 8}) {
          if (is_fixed(j)) continue;
          x_try[j] += norm_dist(rng) * scale * (max_bounds[j] - min_bounds[j]);
        }
      };

      for (int i = 0; i < n_restarts; ++i) {
        Eigen::VectorXd x_try = x0;

        if (i < n_restarts / 3) {
          // Small log-scale perturbations
          for (int j = 0; j < n_log; ++j)
            if (!is_fixed(log_idx[j]))
              x_try[log_idx[j]] += std::log(0.8 + 0.4 * unit_dist(rng));
          perturb_position(x_try, 0.1);
        } else if (i < 2 * n_restarts / 3) {
          // Larger uniform log scale
          const double scale = 0.5 + 1.5 * unit_dist(rng);
          for (int j = 0; j < n_log; ++j)
            if (!is_fixed(log_idx[j])) x_try[log_idx[j]] += std::log(scale);
          perturb_position(x_try, 0.3);
        } else {
          // Random within bounds
          for (int j = 0; j < functor->num_params; ++j)
            if (!is_fixed(j))
              x_try[j] =
                min_bounds[j] + unit_dist(rng) * (max_bounds[j] - min_bounds[j]);
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
          std::cerr << "profile_3d_ic fitting failure: invalid error (" << error
                    << ")\n";
        return false;
      }
      if (I_prf < 1e-7) {
        if (show_error)
          std::cerr << "profile_3d_ic fitting failure: intensity too small (" << I_prf
                    << ")\n";
        return false;
      }
      const double SigX = get_SigX(), SigY = get_SigY(), SigP = get_SigP();
      if (SigX <= 0.0 || SigY <= 0.0) {
        if (show_error)
          std::cerr << "profile_3d_ic fitting failure: non-positive sigma\n";
        return false;
      }
      if (std::abs(SigP) >= 1.0) {
        if (show_error)
          std::cerr << "profile_3d_ic fitting failure: |SigP|>=1 (" << SigP << ")\n";
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
                std::cerr << "profile_3d_ic fitting failure: non-finite model value\n";
              return false;
            }
            if (m > limit) {
              if (show_error)
                std::cerr
                  << "profile_3d_ic fitting failure: model value exceeds limit\n";
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

      // correlation
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
          std::cerr << "profile_3d_ic fitting failure: low correlation (" << corr
                    << ")\n";
        return false;
      }
      if (data_peak > 0.0 && std::abs(profile_peak - data_peak) > 0.25 * data_peak) {
        if (show_error)
          std::cerr << "profile_3d_ic fitting failure: peak height mismatch (profile="
                    << profile_peak << ", data=" << data_peak << ")\n";
        return false;
      }
      return true;
    }
  };

  bool fit_profile_3d_ic(
    scitbx::af::const_ref<vec3<double>, af::c_grid<3>> coords,
    const scitbx::af::versa<double, af::c_grid<3>> intensities,
    const scitbx::af::versa<double, af::c_grid<3>> background_variances,
    TOFProfile3DICParams& profile_params,
    double& I_prf_out,
    boost::optional<scitbx::af::versa<double, af::c_grid<3>>> profile_3d_out =
      boost::none,
    bool update_params = false) {
    const std::array<double, 2> A_bounds = {profile_params.A_min, profile_params.A_max};
    const std::array<double, 2> B_bounds = {profile_params.B_min, profile_params.B_max};
    const std::array<double, 2> R_bounds = {profile_params.R_min, profile_params.R_max};
    const std::array<double, 2> SigX_bounds = {profile_params.SigX_min,
                                               profile_params.SigX_max};
    const std::array<double, 2> SigY_bounds = {profile_params.SigY_min,
                                               profile_params.SigY_max};
    const std::array<double, 2> SigP_bounds = {profile_params.SigP_min,
                                               profile_params.SigP_max};

    TOFProfile3DIC profile(coords,
                           intensities,
                           background_variances,
                           profile_params.A,
                           profile_params.B,
                           profile_params.R,
                           profile_params.SigP,
                           profile_params.HatWidth,
                           profile_params.KConv,
                           A_bounds,
                           B_bounds,
                           R_bounds,
                           SigX_bounds,
                           SigY_bounds,
                           SigP_bounds,
                           profile_params.n_restarts,
                           profile_params.optimize_convolution_params,
                           profile_params.optimize_moderator_params,
                           profile_params.max_drift_factor);

    bool profile_success = true;
    if (profile_params.optimize_profile)
      profile_success = profile.fit(profile_params.show_profile_failures);

    if (profile_success) {
      I_prf_out = profile.calc_intensity();

      if (update_params) {
        profile_params.A = profile.get_A();
        profile_params.B = profile.get_B();
        profile_params.R = profile.get_R();
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

#endif  // DIALS_ALGORITHMS_INTEGRATION_TOF_PROFILE_3D_IC_H

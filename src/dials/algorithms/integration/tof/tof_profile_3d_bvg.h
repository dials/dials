#ifndef DIALS_ALGORITHMS_INTEGRATION_TOF_PROFILE_3D_BVG_H
#define DIALS_ALGORITHMS_INTEGRATION_TOF_PROFILE_3D_BVG_H

#include <array>
#include <algorithm>
#include <cmath>
#include <complex>
#include <scitbx/constants.h>
#include <scitbx/vec3.h>
#include <dials/array_family/scitbx_shared_and_versa.h>

/*
Bivariate Gaussian (BVG) spatial model, and supporting geometry/initial-
estimate helpers, shared by the separable 3D profile fitters (ToF shape ×
BVG(x,y)), e.g. tof_profile_3d_ic.h and tof_profile_3d_ibix.h.
*/

namespace dials { namespace algorithms {

  inline double as_real(double x) {
    return x;
  }
  inline double as_real(const std::complex<double>& x) {
    return x.real();
  }

  // Overflow-guarded exp
  template <typename T>
  inline T exp_safe_t(const T& x) {
    const double r = as_real(x);
    if (r > 700.0) return std::exp(T(700.0));
    if (r < -700.0) return std::exp(T(-700.0));
    return std::exp(x);
  }

  // Bivariate Gaussian PDF at (dx, dy) centred at origin
  template <typename T>
  static T bvg_func_t(double dx, double dy, T SigX, T SigY, T SigP) {
    const double PI = scitbx::constants::pi;
    T rho2 = SigP * SigP;
    if (as_real(rho2) >= 1.0 - 1e-10) rho2 = T(1.0 - 1e-10);
    const T one_minus_rho2 = T(1.0) - rho2;
    const T denom = T(2.0 * PI) * SigX * SigY * std::sqrt(one_minus_rho2);
    if (as_real(denom) < 1e-30) return T(0.0);
    const T z = T(dx * dx) / (SigX * SigX) + T(dy * dy) / (SigY * SigY)
                - T(2.0) * SigP * T(dx * dy) / (SigX * SigY);
    const T val = exp_safe_t(-z / (T(2.0) * one_minus_rho2)) / denom;
    return std::isfinite(as_real(val)) ? val : T(0.0);
  }

  static double bvg_func(double dx, double dy, double SigX, double SigY, double SigP) {
    return bvg_func_t<double>(dx, dy, SigX, SigY, SigP);
  }

  // Recentres shoebox coordinates so the voxel with maximum intensity sits at
  // the origin
  inline scitbx::af::versa<vec3<double>, af::c_grid<3>> get_rel_coords_3d(
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

  // Initial BVG spatial params (from intensity-weighted spatial variance) and
  // per-reflection dXdt/dYdt drift bounds (capping how far the BVG centre may
  // move over the shoebox's ToF range, relative to its spatial extent).
  struct BVGSpatialInit {
    double SigX_init, SigY_init, SigP_init;
    double dXdt_max, dYdt_max;
  };

  inline BVGSpatialInit estimate_bvg_spatial_init(
    const scitbx::af::versa<vec3<double>, af::c_grid<3>>& coords,
    const scitbx::af::versa<double, af::c_grid<3>>& y_norm,
    double SigP,
    const std::array<double, 2>& SigX_bounds,
    const std::array<double, 2>& SigY_bounds,
    const std::array<double, 2>& SigP_bounds,
    double tof_range,
    double max_drift_factor) {
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

    BVGSpatialInit out;
    out.SigX_init = std::min(std::max(SigX_est, SigX_bounds[0]), SigX_bounds[1]);
    out.SigY_init = std::min(std::max(SigY_est, SigY_bounds[0]), SigY_bounds[1]);
    out.SigP_init = std::min(std::max(SigP, SigP_bounds[0]), SigP_bounds[1]);

    const std::size_t nx = coords.accessor()[0];
    const std::size_t ny = coords.accessor()[1];
    const double x_extent = std::abs(coords(nx - 1, 0, 0)[0] - coords(0, 0, 0)[0]);
    const double y_extent = std::abs(coords(0, ny - 1, 0)[1] - coords(0, 0, 0)[1]);
    const double drift_denom = std::max(tof_range, 1e-6);
    out.dXdt_max = max_drift_factor * std::max(x_extent, SigX_bounds[1]) / drift_denom;
    out.dYdt_max = max_drift_factor * std::max(y_extent, SigY_bounds[1]) / drift_denom;
    return out;
  }

}}  // namespace dials::algorithms

#endif  // DIALS_ALGORITHMS_INTEGRATION_TOF_PROFILE_3D_BVG_H

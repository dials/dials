
#ifndef DIALS_ALGORITHMS_SPOT_PREDICTION_SCAN_VARYING_RAY_PREDICTOR_H
#define DIALS_ALGORITHMS_SPOT_PREDICTION_SCAN_VARYING_RAY_PREDICTOR_H

#include <dials/algorithms/spot_prediction/ray_predictor.h>
#include <dials/algorithms/spot_prediction/scan_varying_helpers.h>

namespace dials { namespace algorithms {

  class ScanVaryingRayPredictor {
  public:

    ScanVaryingRayPredictor(
          vec3<double> s0, vec3<double> m2,
          vec2<double> dphi, double dmin)
      : s0_(s0),
        m2_(m2),
        dphi_(dphi),
        s0_mag_(s0.length()),
        dmin_(dmin) {
      DIALS_ASSERT(dmin_ > 0.0);
      dstarmax_ = 1.0 / dmin_;
      dstarmax_sq_ = dstarmax_ * dstarmax_;
    }

    boost::optional<Ray> operator()(const miller_index &h,
        const mat3<double> &A1, const mat3<double> A2, int image, std::size_t step) const {

      // Calculate the reciprocal space vectors
      vec3<double> r1 = A1 * h;
      vec3<double> r2 = A2 * h;
      vec3<double> dr = r2 - r1;
      vec3<double> s0pr1 = s0_ + r1;
      vec3<double> s0pr2 = s0_ + r2;

      // Calculate the distances from the Ewald sphere along radii
      double r1_from_es = s0pr1.length() - s0_mag_;
      double r2_from_es = s0pr2.length() - s0_mag_;

      // Check that the reflection cross the ewald sphere and is within
      // the resolution limit
      bool starts_outside = r1_from_es >= 0.0;
      bool ends_outside = r2_from_es >= 0.0;
      bool is_outside_res_limit = r1.length_sq() > dstarmax_sq_;
      if (starts_outside == ends_outside || is_outside_res_limit) {
        return boost::optional<Ray>();
      }

      // Solve the equation |s0 + r1 + alpha * dr| = |s0| for alpha. This is
      // equivalent to solving the quadratic equation
      //
      // alpha^2*dr.dr + 2*alpha(s0 + r1).dr + 2*s0.r1 + r1.r1 = 0
      af::small<double, 2> roots = reeke_detail::solve_quad(
          dr.length_sq(),
          2.0 * s0pr1 * dr,
          r1.length_sq() + 2.0*s0_*r1);

      // Choose a root that lies in [0,1]
      double alpha;
      if (0.0 <= roots[0] && roots[0] <= 1.0) {
        alpha = roots[0];
      } else if (0.0 <= roots[1] && roots[1] <= 1.0) {
        alpha = roots[1];
      } else {
        DIALS_ASSERT(false);
      }

      // Calculate the scattering vector and rotation angle
      vec3<double> s1 = r1 + alpha * dr + s0_;
      double angle = dphi_[0] + (image + alpha * step) * dphi_[1];

      // Return the ray
      return Ray(s1, angle, starts_outside);
    }


  private:
    vec3<double> s0_;
    vec3<double> m2_;
    vec2<double> dphi_;
    double s0_mag_;
    double dmin_;
    double dstarmax_;
    double dstarmax_sq_;
  };

}} // namespace dials::algorithms

#endif // DIALS_ALGORITHMS_SPOT_PREDICTION_SCAN_VARYING_RAY_PREDICTOR_H

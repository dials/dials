/*
 * scan_varying_ray_predictor.h
 *
 *  Copyright (C) 2013 Diamond Light Source, CCP4
 *
 *  Author: David Waterman
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */

#ifndef DIALS_ALGORITHMS_SPOT_PREDICTION_SCAN_VARYING_RAY_PREDICTOR_H
#define DIALS_ALGORITHMS_SPOT_PREDICTION_SCAN_VARYING_RAY_PREDICTOR_H

#include <scitbx/vec2.h>
#include <scitbx/vec3.h>
#include <scitbx/mat3.h>
#include <cctbx/miller.h>
#include <dials/model/data/ray.h>
#include <dials/algorithms/spot_prediction/scan_varying_helpers.h>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/error.h>

namespace dials { namespace algorithms {

  using dials::model::Ray;
  using scitbx::mat3;
  using scitbx::vec2;
  using scitbx::vec3;

  /**
   * Reflection prediction for a relp within a small segment of a scan, which we
   * assume to be a single image.
   *
   * The path of the relp through reciprocal space between the start and end of
   * the image is approximated by a general linear transformation (not just a
   * rotation).
   *
   * Temporarily, we initialise with a list of N+1 UB matrices, where N is the
   * total number of images. In future we will simple pass a crystal model,
   * which will store its own per-image UB matrix.
   *
   * The first operator() method assumes only the crystal model varies with
   * image number, whilst the other models remain static. The overload allows
   * also the beam to differ between the start and end of the image.
   */
  class ScanVaryingRayPredictor {
  public:
    // Typedef the miller_index type
    typedef cctbx::miller::index<> miller_index;

    /**
     * Initialise the predictor.
     * @param s0 The beam vector
     * @param m2 The rotation axis
     * @param dphi The oscillation (phi0, dphi)
     * @param dmin The resolution limit
     */
    ScanVaryingRayPredictor(vec3<double> s0,
                            vec3<double> m2,
                            int frame0,
                            vec2<double> dphi,
                            double dmin)
        : s0_(s0),
          m2_(m2.normalize()),
          frame0_(frame0),
          dphi_(dphi),
          s0_mag_(s0.length()),
          dmin_(dmin) {
      DIALS_ASSERT(std::abs(dphi_[1]) > 0.0);
      DIALS_ASSERT(s0_mag_ > 0.0);
      DIALS_ASSERT(dmin_ > 0.0);
      dstarmax_ = 1.0 / dmin_;
      dstarmax_sq_ = dstarmax_ * dstarmax_;
    }

    /**
     * Predict the ray for the given Miller index on the given image, where the
     * UB matrix differs between the start and end of the step
     * @param h The miller index
     * @param A1 The setting matrix for the beginning of the step.
     * @param A2 The setting matrix for the end of the step.
     * @param image The image index
     * @param step The step to predict over.
     * @returns The ray if predicted
     */
    boost::optional<Ray> operator()(const miller_index &h,
                                    const mat3<double> &A1,
                                    const mat3<double> &A2,
                                    int image,
                                    std::size_t step) const {
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
        dr.length_sq(), 2.0 * s0pr1 * dr, r1.length_sq() + 2.0 * s0_ * r1);

      // Choose a root that lies in [0,1]
      double alpha;
      if (0.0 <= roots[0] && roots[0] <= 1.0) {
        alpha = roots[0];
      } else if (0.0 <= roots[1] && roots[1] <= 1.0) {
        alpha = roots[1];
      } else {
        return boost::optional<Ray>();
      }

      // Calculate the scattering vector and rotation angle
      vec3<double> s1 = r1 + alpha * dr + s0_;
      double angle = dphi_[0] + (image - frame0_ + alpha * step) * dphi_[1];

      // Return the ray
      return Ray(s1, angle, starts_outside);
    }

    /**
     * Predict the ray for the given Miller index on the given image, where the
     * UB matrix and the s0 vector differs between the start and end of the
     * step.
     * @param h The miller index
     * @param A1 The setting matrix for the beginning of the step.
     * @param A2 The setting matrix for the end of the step.
     * @param s0a The s0 vector for the beginning of the step.
     * @param s0b The s0 vector for the end of the step.
     * @param image The image index
     * @param step The step to predict over.
     * @returns The ray if predicted
     */
    boost::optional<Ray> operator()(const miller_index &h,
                                    const mat3<double> &A1,
                                    const mat3<double> &A2,
                                    const vec3<double> &s0a,
                                    const vec3<double> &s0b,
                                    int image,
                                    std::size_t step) const {
      // Calculate the reciprocal space vectors
      vec3<double> r1 = A1 * h;
      vec3<double> r2 = A2 * h;
      vec3<double> dr = r2 - r1;
      vec3<double> s0pr1 = s0a + r1;
      vec3<double> s0pr2 = s0b + r2;

      // Calculate the distances from the Ewald spheres along radii
      double s0a_mag = s0a.length();
      double s0b_mag = s0b.length();
      double r1_from_es = s0pr1.length() - s0a_mag;
      double r2_from_es = s0pr2.length() - s0b_mag;

      // Check that the reflection cross the ewald sphere and is within
      // the resolution limit
      bool starts_outside = r1_from_es >= 0.0;
      bool ends_outside = r2_from_es >= 0.0;
      bool is_outside_res_limit = r1.length_sq() > dstarmax_sq_;
      if (starts_outside == ends_outside || is_outside_res_limit) {
        return boost::optional<Ray>();
      }

      // Calculate distance of r1 from the start Ewald sphere along the
      // direction of linear change dr. That implies solving the equation
      // |s0a + r1 + alpha * dr| = |s0a| for alpha. This is equivalent to
      // solving the quadratic equation
      //
      // alpha^2*dr.dr + 2*alpha(s0a + r1).dr + 2*s0a.r1 + r1.r1 = 0
      af::small<double, 2> roots = reeke_detail::solve_quad(
        dr.length_sq(), 2.0 * s0pr1 * dr, r1.length_sq() + 2.0 * s0a * r1);

      // Choose a root that lies in [0,1]
      double alpha1;
      if (0.0 <= roots[0] && roots[0] <= 1.0) {
        alpha1 = roots[0];
      } else if (0.0 <= roots[1] && roots[1] <= 1.0) {
        alpha1 = roots[1];
      } else {
        return boost::optional<Ray>();
      }

      // Now calculate the distance of r2 from the end Ewald sphere along the
      // direction of linear change -1.0*dr. That implies solving the equation
      // |s0b + r2 - alpha * dr| = |s0b| for alpha. This is equivalent to
      // solving the quadratic equation
      //
      // alpha^2*dr.dr - 2*alpha(s0b + r2).dr + 2*s0b.r2 + r2.r2 = 0
      roots = reeke_detail::solve_quad(
        dr.length_sq(), -2.0 * s0pr2 * dr, r2.length_sq() + 2.0 * s0b * r2);
      // Choose a root that lies in [0,1]
      double alpha2;
      if (0.0 <= roots[0] && roots[0] <= 1.0) {
        alpha2 = roots[0];
      } else if (0.0 <= roots[1] && roots[1] <= 1.0) {
        alpha2 = roots[1];
      } else {
        return boost::optional<Ray>();
      }

      // Calculate alpha, the fraction along the linear step, as the distance
      // from the Ewald sphere at the start compared to the total distance
      // travelled relative to the Ewald sphere
      double alpha = alpha1 / (alpha1 + alpha2);

      // Linear approximation to the s0 vector at intersection
      vec3<double> us0a = s0a.normalize();
      vec3<double> us0_at_intersection = alpha * (s0b.normalize() - us0a) + us0a;
      double wavenumber = (s0a_mag + s0b_mag) * 0.5;
      vec3<double> s0_at_intersection = wavenumber * us0_at_intersection;

      // Calculate the scattering vector and rotation angle
      vec3<double> s1 = r1 + alpha * dr + s0_at_intersection;
      double angle = dphi_[0] + (image - frame0_ + alpha * step) * dphi_[1];

      // Return the ray
      return Ray(s1, angle, starts_outside);
    }

  private:
    vec3<double> s0_;
    vec3<double> m2_;
    int frame0_;
    vec2<double> dphi_;
    double s0_mag_;
    double dmin_;
    double dstarmax_;
    double dstarmax_sq_;
  };

}}  // namespace dials::algorithms

#endif  // DIALS_ALGORITHMS_SPOT_PREDICTION_SCAN_VARYING_RAY_PREDICTOR_H

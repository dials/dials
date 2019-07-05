/*
 * ray_predictor.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_SPOT_PREDICTION_RAY_PREDICTOR_H
#define DIALS_ALGORITHMS_SPOT_PREDICTION_RAY_PREDICTOR_H

#include <scitbx/constants.h>
#include <scitbx/vec2.h>
#include <scitbx/vec3.h>
#include <scitbx/mat3.h>
#include <cctbx/miller.h>
#include <dxtbx/model/scan_helpers.h>
#include <dials/model/data/ray.h>
#include "rotation_angles.h"

namespace dials { namespace algorithms {

  // Using lots of stuff from other namespaces
  using dxtbx::model::is_angle_in_range;
  using dxtbx::model::mod_2pi;
  using model::Ray;
  using scitbx::mat3;
  using scitbx::rad_as_deg;
  using scitbx::vec2;
  using scitbx::vec3;

  class ScanStaticRayPredictor {
  public:
    // Typedef the miller_index type
    typedef cctbx::miller::index<> miller_index;

    /**
     * Initialise the ray predictor.
     * @param s0 The incident beam vector
     * @param m2 The rotation axis
     * @param dphi The total oscillation range
     */
    ScanStaticRayPredictor(vec3<double> s0,
                           vec3<double> m2,
                           mat3<double> fixed_rotation,
                           mat3<double> setting_rotation,
                           vec2<double> dphi)
        : calculate_rotation_angles_(setting_rotation.inverse() * s0, m2),
          fixed_rotation_(fixed_rotation),
          setting_rotation_(setting_rotation),
          dphi_(dphi),
          s0_(s0),
          m2_(m2.normalize()),
          s0_m2_plane(s0.cross(setting_rotation * m2).normalize()) {}

    /** Virtual destructor to allow inheritance */
    virtual ~ScanStaticRayPredictor() {}

    /**
     * Predict the spot locations on the image detector.
     *
     * The algorithm performs the following procedure:
     *
     *  - For the miller index, the rotation angle at which the diffraction
     *    conditions are met is calculated.
     *
     *  - The rotation angles are then checked to see if they are within the
     *    rotation range.
     *
     *  - The reciprocal lattice vectors are then calculated, followed by the
     *    diffracted beam vector for each reflection.
     *
     *  - The diffracted beam vectors are then classified according to whether
     *    the reciprocal lattice point is entering or exiting the Ewald sphere.
     *
     * @param h The miller index
     * @returns An array of predicted reflections
     */
    af::small<Ray, 2> operator()(miller_index h, mat3<double> UB) const {
      af::small<Ray, 2> rays;

      // Calculate the reciprocal space vector
      vec3<double> pstar0 = fixed_rotation_ * UB * h;

      // Try to calculate the diffracting rotation angles
      vec2<double> phi;
      try {
        phi = calculate_rotation_angles_(pstar0);
      } catch (error) {
        return rays;
      }

      // Loop through the 2 rotation angles
      for (std::size_t i = 0; i < phi.size(); ++i) {
        // Check that the angles are within the rotation range
        if (!is_angle_in_range(dphi_, phi[i])) {
          continue;
        }

        // Calculate the reciprocal space vector and diffracted beam vector
        vec3<double> pstar =
          setting_rotation_ * pstar0.unit_rotate_around_origin(m2_, phi[i]);
        vec3<double> s1 = s0_ + pstar;

        double small = 1.0e-8;

        DIALS_ASSERT(std::abs(s1.length() - s0_.length()) < small);

        // Calculate the direction of reflection passage
        bool entering = s1 * s0_m2_plane < 0.;

        // Add the reflection
        rays.push_back(Ray(s1, mod_2pi(phi[i]), entering));
      }
      return rays;
    }

    af::small<Ray, 2> from_reciprocal_lattice_vector(vec3<double> pstar0) const {
      af::small<Ray, 2> rays;

      // Try to calculate the diffracting rotation angles
      vec2<double> phi;
      try {
        phi = calculate_rotation_angles_(pstar0);
      } catch (error) {
        return rays;
      }

      // Loop through the 2 rotation angles
      for (std::size_t i = 0; i < phi.size(); ++i) {
        // Check that the angles are within the rotation range
        if (!is_angle_in_range(dphi_, phi[i])) {
          continue;
        }

        // Calculate the reciprocal space vector and diffracted beam vector
        vec3<double> pstar =
          setting_rotation_ * pstar0.unit_rotate_around_origin(m2_, phi[i]);
        vec3<double> s1 = s0_ + pstar;

        double small = 1.0e-8;

        DIALS_ASSERT(std::abs(s1.length() - s0_.length()) < small);

        // Calculate the direction of reflection passage
        bool entering = s1 * s0_m2_plane < 0.;

        // Add the reflection
        rays.push_back(Ray(s1, mod_2pi(phi[i]), entering));
      }
      return rays;
    }

  private:
    RotationAngles calculate_rotation_angles_;
    mat3<double> fixed_rotation_;
    mat3<double> setting_rotation_;
    vec2<double> dphi_;
    vec3<double> s0_;
    vec3<double> m2_;
    vec3<double> s0_m2_plane;
  };

  // Typedef the miller_index type
  typedef cctbx::miller::index<> miller_index;

  class RayPredictor {
  public:
    // Typedef the miller_index type
    typedef cctbx::miller::index<> miller_index;

    /**
     * Initialise the ray predictor.
     * @param s0 The incident beam vector
     * @param m2 The rotation axis
     * @param dphi The total oscillation range
     */
    RayPredictor(vec3<double> m2,
                 mat3<double> fixed_rotation,
                 mat3<double> setting_rotation)
        : fixed_rotation_(fixed_rotation),
          setting_rotation_(setting_rotation),
          setting_rotation_inv_(setting_rotation_.inverse()),
          m2_(m2.normalize()) {}

    /** Virtual destructor to allow inheritance */
    virtual ~RayPredictor() {}

    /**
     * Predict the spot locations on the image detector.
     *
     * The algorithm performs the following procedure:
     *
     *  - For the miller index, the rotation angle at which the diffraction
     *    conditions are met is calculated.
     *
     *  - The rotation angles are then checked to see if they are within the
     *    rotation range.
     *
     *  - The reciprocal lattice vectors are then calculated, followed by the
     *    diffracted beam vector for each reflection.
     *
     *  - The diffracted beam vectors are then classified according to whether
     *    the reciprocal lattice point is entering or exiting the Ewald sphere.
     *
     * @param h The miller index
     * @returns An array of predicted reflections
     */
    Ray operator()(mat3<double> UB,
                   vec3<double> s0,
                   vec3<double> h,
                   bool entering_in) const {
      vec3<double> s0_m2_plane = s0.cross(setting_rotation_ * m2_).normalize();

      // Calculate the reciprocal space vector
      vec3<double> pstar0 = fixed_rotation_ * UB * h;

      // Try to calculate the diffracting rotation angles
      vec2<double> phi = rotation_angles(setting_rotation_inv_ * s0, m2_, pstar0);

      // Loop through the 2 rotation angles
      for (std::size_t i = 0; i < phi.size(); ++i) {
        // Calculate the reciprocal space vector and diffracted beam vector
        vec3<double> pstar =
          setting_rotation_ * pstar0.unit_rotate_around_origin(m2_, phi[i]);
        vec3<double> s1 = s0_ + pstar;

        double small = 1.0e-8;

        DIALS_ASSERT(std::abs(s1.length() - s0_.length()) < small);

        // Calculate the direction of reflection passage
        bool entering = s1 * s0_m2_plane < 0.;
        if (entering == entering_in) {
          return Ray(s1, mod_2pi(phi[i]), entering);
        }
      }
      throw DIALS_ERROR("No ray");
      return Ray();
    }

  private:
    mat3<double> fixed_rotation_;
    mat3<double> setting_rotation_;
    mat3<double> setting_rotation_inv_;
    vec3<double> s0_;
    vec3<double> m2_;
    vec3<double> s0_m2_plane;
  };

}}  // namespace dials::algorithms

#endif  // DIALS_ALGORITHMS_SPOT_PREDICTION_RAY_PREDICTOR_H

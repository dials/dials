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
#include <dials/model/data/reflection.h>
#include "rotation_angles.h"

namespace dials { namespace algorithms {

  // Using lots of stuff from other namespaces
  using scitbx::rad_as_deg;
  using scitbx::vec2;
  using scitbx::vec3;
  using scitbx::mat3;
  using dxtbx::model::mod_2pi;
  using dxtbx::model::is_angle_in_range;
  using model::Reflection;

  struct Ray {
    vec3<double> s1;
    double angle;
    bool entering;

    Ray() {}

    Ray(vec3<double> s1_, double angle_, bool entering_)
      : s1(s1_),
        angle(angle_),
        entering(entering_) {}
  };

  class RayPredictor2 {
  public:

    // Typedef the miller_index type
    typedef cctbx::miller::index <> miller_index;

    /**
     * Initialise the ray predictor.
     * @param s0 The incident beam vector
     * @param m2 The rotation axis
     * @param dphi The total oscillation range
     */
    RayPredictor2(vec3 <double> s0, vec3 <double> m2,
                 vec2 <double> dphi)
      : calculate_rotation_angles_(s0, m2),
        dphi_(dphi),
        s0_(s0),
        m2_(m2.normalize()),
        s0_m2_plane(s0.cross(m2).normalize()){}

    /** Virtual destructor to allow inheritance */
    virtual ~RayPredictor2() {}

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
    af::small<Ray,2> operator()(miller_index h, mat3 <double> UB) const {

      af::small<Ray,2> rays;

      // Calculate the reciprocal space vector
      vec3 <double> pstar0 = UB * h;

      // Try to calculate the diffracting rotation angles
      vec2 <double> phi;
      try {
        phi = calculate_rotation_angles_(pstar0);
      } catch(error) {
        return rays;
      }

      // Loop through the 2 rotation angles
      for (std::size_t i = 0; i < phi.size(); ++i) {

        // Check that the angles are within the rotation range
        if (!is_angle_in_range(dphi_, phi[i])) {
          continue;
        }

        // Calculate the reciprocal space vector and diffracted beam vector
        vec3 <double> pstar = pstar0.unit_rotate_around_origin(m2_, phi[i]);
        vec3 <double> s1 = s0_ + pstar;

        // Calculate the direction of reflection passage
        bool entering = s1 * s0_m2_plane < 0.;

        // Add the reflection
        rays.push_back(Ray(s1, mod_2pi(phi[i]), entering));
      }
      return rays;
    }

  private:
    RotationAngles calculate_rotation_angles_;
    vec2 <double> dphi_;
    vec3 <double> s0_;
    vec3 <double> m2_;
    vec3 <double> s0_m2_plane;
  };


  // Typedef the miller_index type
  typedef cctbx::miller::index <> miller_index;

  /** A class to perform spot prediction. */
  class RayPredictor {
  public:

    // A load of useful typedefs
    typedef Reflection reflection_type;
    typedef af::shared <reflection_type> reflection_list_type;

    /**
     * Initialise the ray predictor.
     * @param s0 The incident beam vector
     * @param m2 The rotation axis
     * @param dphi The total oscillation range
     */
    RayPredictor(vec3 <double> s0, vec3 <double> m2,
                 vec2 <double> dphi)
      : calculate_rotation_angles_(s0, m2),
        dphi_(dphi),
        s0_(s0),
        m2_(m2.normalize()),
        s0_m2_plane(s0.cross(m2).normalize()){}

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
    reflection_list_type
    operator()(miller_index h, mat3 <double> UB) const {

      reflection_list_type reflections;

      // Calculate the reciprocal space vector
      vec3 <double> pstar0 = UB * h;

      // Try to calculate the diffracting rotation angles
      vec2 <double> phi;
      try {
        phi = calculate_rotation_angles_(pstar0);
      } catch(error) {
        return reflections;
      }

      // Loop through the 2 rotation angles
      for (std::size_t i = 0; i < phi.size(); ++i) {

        // Check that the angles are within the rotation range
        if (!is_angle_in_range(dphi_, phi[i])) {
          continue;
        }

        // Calculate the reciprocal space vector and diffracted beam vector
        vec3 <double> pstar = pstar0.unit_rotate_around_origin(m2_, phi[i]);
        vec3 <double> s1 = s0_ + pstar;

        // Calculate the direction of reflection passage
        bool entering = s1 * s0_m2_plane < 0.;

        // Add the reflection
        reflections.push_back(reflection_type(h, mod_2pi(phi[i]), s1, entering));
      }
      return reflections;
    }

    /**
     * For a given set of miller indices and a single UB matrix, predict
     * the detector coordinates.
     * @param miller_indices The array of miller indices.
     */
    reflection_list_type
    operator()(const af::const_ref<miller_index> &miller_indices,
               mat3 <double> UB) const {
      reflection_list_type reflections;
      for (std::size_t i = 0; i < miller_indices.size(); ++i) {
        reflection_list_type r = operator()(miller_indices[i], UB);
        for (std::size_t j = 0; j < r.size(); ++j) {
          reflections.push_back(r[j]);
        }
      }
      return reflections;
    }

  private:

    RotationAngles calculate_rotation_angles_;
    vec2 <double> dphi_;
    vec3 <double> s0_;
    vec3 <double> m2_;
    vec3 <double> s0_m2_plane;
  };

}} // namespace dials::algorithms

#endif // DIALS_ALGORITHMS_SPOT_PREDICTION_RAY_PREDICTOR_H

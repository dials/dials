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
#ifndef DIALS_ALGORITHMS_RAY_PREDICTION_RAY_PREDICTOR_H
#define DIALS_ALGORITHMS_RAY_PREDICTION_RAY_PREDICTOR_H

#include <scitbx/constants.h>
#include <scitbx/array_family/small.h>
#include <scitbx/array_family/flex_types.h>
#include <dials/model/experiment/scan_helpers.h>
#include "index_generator.h"
#include "rotation_angles.h"

namespace dials { namespace algorithms {

  // Using lots of stuff from other namespaces
  using scitbx::rad_as_deg;
  using scitbx::vec2;
  using scitbx::vec3;
  using scitbx::mat3;
  using scitbx::af::small;
  using scitbx::af::flex_double;
  using model::is_scan_angle_valid;
  using model::mod_360;

  // Typedef the miller_index and flex_miller_index types
  typedef cctbx::miller::index <> miller_index;
  typedef scitbx::af::flex <miller_index> ::type flex_miller_index;

  /** A class to perform spot prediction. */
  template <typename BeamType,
            typename GoniometerType,
            typename ScanType,
            typename ReflectionType>
  class RayPredictor {
  public:

    // A load of useful typedefs
    typedef BeamType beam_type;
    typedef GoniometerType goniometer_type;
    typedef ScanType scan_type;
    typedef ReflectionType reflection_type;
    typedef scitbx::af::shared <reflection_type> reflection_list_type;

    /**
     * Initialise the spot predictor.
     * @param beam The beam parameters
     * @param gonio The goniometer parameters
     * @param scan The scan parameters
     * @param unit_cell The unit cell parameters
     * @param space_group_type The space group struct
     * @param UB The ub matrix
     * @param d_min The resolution
     */
    RayPredictor(const beam_type &beam,
                 const goniometer_type &gonio,
                 const scan_type &scan,
                 mat3 <double> UB)
      : calculate_rotation_angles_(
          beam.get_direction(),
          gonio.get_rotation_axis()),
        is_angle_valid_(scan),
        UB_(gonio.get_fixed_rotation() * UB),
        s0_(beam.get_direction()),
        m2_(gonio.get_rotation_axis()) {}

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
     * @param h The miller index
     * @returns An array of predicted reflections
     */
    reflection_list_type
    operator()(miller_index h) const {

      reflection_list_type reflections;

      // Calculate the reciprocal space vector
      vec3 <double> pstar0 = UB_ * h;

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
        double phi_deg = mod_360(rad_as_deg(phi[i]));
        if (!is_angle_valid_(phi_deg)) {
          continue;
        }

        // Calculate the reciprocal space vector and diffracted beam vector
        vec3 <double> pstar = pstar0.unit_rotate_around_origin(m2_, phi[i]);
        vec3 <double> s1 = s0_ + pstar;
      
        // Add the reflection
        reflections.push_back(reflection_type(h, phi_deg, s1));
      }
      return reflections;
    }

    /**
     * For a given set of miller indices, predict the detector coordinates.
     * @param miller_indices The array of miller indices.
     */
    reflection_list_type
    operator()(const flex_miller_index &miller_indices) const {
      reflection_list_type reflections;
      for (std::size_t i = 0; i < miller_indices.size(); ++i) {
        reflection_list_type r = operator()(miller_indices[i]);
        for (std::size_t j = 0; j < r.size(); ++j) {
          reflections.push_back(r[j]);
        }
      }
      return reflections;
    }

  private:

    RotationAngles calculate_rotation_angles_;
    is_scan_angle_valid <scan_type> is_angle_valid_;
    mat3 <double> UB_;
    vec3 <double> s0_;
    vec3 <double> m2_;
  };

}} // namespace dials::algorithms

#endif // DIALS_ALGORITHMS_RAY_PREDICTION_RAY_PREDICTOR_H

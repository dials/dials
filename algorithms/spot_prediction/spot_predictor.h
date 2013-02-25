/*
 * spot_predictor.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_SPOT_PREDICTION_SPOT_PREDICTOR_H
#define DIALS_ALGORITHMS_SPOT_PREDICTION_SPOT_PREDICTOR_H

#include <scitbx/constants.h>
#include <scitbx/array_family/flex_types.h>
#include <dials/model/experiment/scan_helpers.h>
#include <dials/model/experiment/detector_helpers.h>
#include "index_generator.h"
#include "rotation_angles.h"

namespace dials { namespace algorithms {

  // Using lots of stuff from other namespaces
  using scitbx::rad_as_deg;
  using scitbx::vec2;
  using scitbx::vec3;
  using scitbx::mat3;
  using scitbx::af::flex_double;
  using model::is_scan_angle_valid;
  using model::diffracted_beam_intersection;
  using model::get_all_frames_from_angle;
  using model::mod_2pi;

  // Typedef the miller_index and flex_miller_index types
  typedef cctbx::miller::index <> miller_index;
  typedef scitbx::af::flex <miller_index> ::type flex_miller_index;

  /** A class to perform spot prediction. */
  template <typename BeamType,
            typename DetectorType,
            typename GoniometerType,
            typename ScanType,
            typename ReflectionType>
  class SpotPredictor {
  public:

    // A load of useful typedefs
    typedef BeamType beam_type;
    typedef DetectorType detector_type;
    typedef GoniometerType goniometer_type;
    typedef ScanType scan_type;
    typedef ReflectionType reflection_type;
    typedef scitbx::af::shared <reflection_type> reflection_list_type;
    typedef typename detector_type::coordinate_type detector_coordinate_type;

    /**
     * Initialise the spot predictor.
     * @param beam The beam parameters
     * @param detector The detector parameters
     * @param gonio The goniometer parameters
     * @param scan The scan parameters
     * @param unit_cell The unit cell parameters
     * @param space_group_type The space group struct
     * @param ub_matrix The ub matrix
     * @param d_min The resolution
     */
    SpotPredictor(const beam_type &beam,
                  const detector_type &detector,
                  const goniometer_type &gonio,
                  const scan_type &scan,
                  const cctbx::uctbx::unit_cell &unit_cell,
                  const cctbx::sgtbx::space_group_type &space_group,
                  mat3 <double> ub_matrix,
                  double d_min)
      : index_generator_(unit_cell, space_group, false, d_min),
        calculate_rotation_angles_(
          beam.get_direction(),
          gonio.get_rotation_axis()),
        is_angle_valid_(scan),
        get_detector_coord_(detector),
        get_frame_numbers_(scan),
        beam_(beam),
        detector_(detector),
        gonio_(gonio),
        scan_(scan),
        ub_matrix_(gonio.get_fixed_rotation() * ub_matrix),
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
     *  - The image volume coordinates are then calculated for each reflection.
     *
     *  - The image volume coordinates are then checked to see if they are
     *    within the image volume itself.
     *
     * @param h The miller index
     * @returns An array of predicted reflections
     */
    reflection_list_type
    operator()(miller_index h) const {

      reflection_list_type reflections;

      // Calculate the reciprocal space vector
      vec3 <double> pstar0 = ub_matrix_ * h;

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
        phi[i] = mod_2pi(phi[i]);
        if (!is_angle_valid_(phi[i])) {
          continue;
        }

        // Calculate the reciprocal space vector and diffracted beam vector
        vec3 <double> pstar = pstar0.unit_rotate_around_origin(m2_, phi[i]);
        vec3 <double> s1 = s0_ + pstar;

        // Try to calculate the detector coordinate
        detector_coordinate_type coord;
        try {
          coord = get_detector_coord_(s1);
        } catch(error) {
          continue;
        }

        // Get the list of frames at which the reflection will be observed
        // and add the predicted observations to the list of reflections
        flex_double frames = get_frame_numbers_(phi[i]);
        for (std::size_t j = 0; j < frames.size(); ++j) {
          reflection_type r = reflection_type(
            h, phi[i], s1, coord, frames[j]);
          reflections.push_back(r);
        }
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

    /**
     * Generate a set of miller indices and predict the detector coordinates.
     */
    reflection_list_type
    operator()() {
      reflection_list_type reflections;

      // Continue looping until we run out of miller indices
      for (;;) {

        // Get the next miller index
        miller_index h = index_generator_.next();
        if (h.is_zero()) {
          break;
        }

        // Predict the spot location for the miller index
        reflection_list_type r = operator()(h);
        for (std::size_t j = 0; j < r.size(); ++j) {
          reflections.push_back(r[j]);
        }
      }
      return reflections;
    }

  private:

    IndexGenerator index_generator_;
    RotationAngles calculate_rotation_angles_;
    is_scan_angle_valid <scan_type> is_angle_valid_;
    diffracted_beam_intersection <detector_type> get_detector_coord_;
    get_all_frames_from_angle <scan_type> get_frame_numbers_;
    beam_type beam_;
    detector_type detector_;
    goniometer_type gonio_;
    scan_type scan_;
    mat3 <double> ub_matrix_;
    vec3 <double> s0_;
    vec3 <double> m2_;
  };

}} // namespace dials::algorithms

#endif // DIALS_ALGORITHMS_SPOT_PREDICTION_SPOT_PREDICTOR_H

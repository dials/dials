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
#include <dials/model/experiment/beam.h>
#include <dials/model/experiment/scan.h>
#include <dials/model/experiment/detector.h>
#include <dials/model/experiment/goniometer.h>
#include <dials/model/experiment/scan_helpers.h>
#include <dials/model/experiment/detector_helpers.h>
#include <dials/model/data/reflection.h>
#include "index_generator.h"
#include "rotation_angles.h"

namespace dials { namespace algorithms {

  // Using lots of stuff from other namespaces
  using scitbx::rad_as_deg;
  using scitbx::vec2;
  using scitbx::vec3;
  using scitbx::mat3;
  using scitbx::af::flex_double;
  using model::Beam;
  using model::Scan;
  using model::FlatPanelDetector;
  using model::Goniometer;
  using model::Reflection;
  using model::ReflectionList;
  using model::is_scan_angle_valid;
  using model::is_coordinate_valid;
  using model::diffracted_beam_detector_coord;
  using model::get_all_frames_from_angle;
  using model::mod_360;

  // Typedef the miller_index and flex_miller_index types
  typedef cctbx::miller::index <> miller_index;
  typedef scitbx::af::flex <miller_index> ::type flex_miller_index;

  /** A class to perform spot prediction. */
  class SpotPredictor {

  public:

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
    SpotPredictor(const Beam &beam,
                  const FlatPanelDetector &detector,
                  const Goniometer &gonio,
                  const Scan &scan,
                  const cctbx::uctbx::unit_cell &unit_cell,
                  const cctbx::sgtbx::space_group_type &space_group_type,
                  mat3 <double> ub_matrix,
                  double d_min)
      : index_generator_(unit_cell, space_group_type, false, d_min),
        calculate_rotation_angles_(
          beam.get_direction(),
          gonio.get_rotation_axis()),
        is_angle_valid_(scan),
        is_coord_valid_(detector),
        get_image_coord_(detector),
        get_frame_numbers_(scan),
        beam_(beam),
        detector_(detector),
        gonio_(gonio),
        scan_(scan),
        ub_matrix_(gonio.get_fixed_rotation() * ub_matrix),
        s0_(beam.get_direction()),
        m2_(gonio.get_rotation_axis().normalize()) {}

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
    scitbx::af::shared <Reflection> predict(miller_index h) const {

      scitbx::af::shared <Reflection> reflections;

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

        // Convert angle to degrees
        double phi_deg = mod_360(rad_as_deg(phi[i]));

        // Check that the angles are within the rotation range
        if (!is_angle_valid_(phi_deg)) {
          continue;
        }

        // Calculate the reciprocal space vector and diffracted beam vector
        vec3 <double> pstar = pstar0.unit_rotate_around_origin(m2_, phi[i]);
        vec3 <double> s1 = s0_ + pstar;

        // Try to calculate the detector coordinate
        vec2 <double> xy;
        try {
          xy = get_image_coord_(s1);
        } catch(error) {
          continue;
        }

        // Check the detector coordinate is valid and add the
        // elements to the arrays. NB. up to now, we have used
        // angles in radians, convert them to degrees before adding
        // them to the rotation angle array.
        if (!is_coord_valid_(xy)) {
          continue;
        }

        // Get the list of frames at which the reflection will be observed
        // and add the predicted observations to the list of reflections
        flex_double frames = get_frame_numbers_(phi[i]);
        for (std::size_t j = 0; j < frames.size(); ++j) {
          reflections.push_back(Reflection(h, phi_deg, s1, xy, frames[j]));
        }
      }
      return reflections;
    }

    /**
     * For a given set of miller indices, predict the detector coordinates.
     * @param miller_indices The array of miller indices.
     */
    scitbx::af::shared <Reflection>
    predict(const flex_miller_index &miller_indices) const {
      scitbx::af::shared <Reflection> reflections;
      for (std::size_t i = 0; i < miller_indices.size(); ++i) {
        scitbx::af::shared <Reflection> r = predict(miller_indices[i]);
        for (std::size_t j = 0; j < r.size(); ++j) {
          reflections.push_back(r[j]);
        }
      }
      return reflections;
    }

    /**
     * Generate a set of miller indices and predict the detector coordinates.
     */
    scitbx::af::shared <Reflection>
    predict() {
      scitbx::af::shared <Reflection> reflections;

      // Continue looping until we run out of miller indices
      for (;;) {

        // Get the next miller index
        miller_index h = index_generator_.next();
        if (h.is_zero()) {
          break;
        }

        // Predict the spot location for the miller index
        scitbx::af::shared <Reflection> r = predict(h);
        for (std::size_t j = 0; j < r.size(); ++j) {
          reflections.push_back(r[j]);
        }
      }
      return reflections;
    }

  private:

      IndexGenerator index_generator_;
      RotationAngles calculate_rotation_angles_;
      is_scan_angle_valid <Scan> is_angle_valid_;
      is_coordinate_valid <FlatPanelDetector> is_coord_valid_;
      diffracted_beam_detector_coord <FlatPanelDetector> get_image_coord_;
      get_all_frames_from_angle <Scan> get_frame_numbers_;
      Beam beam_;
      FlatPanelDetector detector_;
      Goniometer gonio_;
      Scan scan_;
      mat3 <double> ub_matrix_;
      vec3 <double> s0_;
      vec3 <double> m2_;
  };

}} // namespace dials::algorithms

#endif // DIALS_ALGORITHMS_SPOT_PREDICTION_SPOT_PREDICTOR_H

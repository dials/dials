/*
 * multi_panel_spot_predictor.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_ALGORITHMS_SPOT_PREDICTION_MULTI_PANEL_SPOT_PREDICTOR_H
#define DIALS_ALGORITHMS_SPOT_PREDICTION_MULTI_PANEL_SPOT_PREDICTOR_H

#include <scitbx/constants.h>
#include <dxtbx/model/scan_helpers.h>
#include <dxtbx/model/scan.h>
#include <dxtbx/model/beam.h>
#include <dxtbx/model/goniometer.h>
#include <dxtbx/model/multi_panel_detector.h>
#include <dials/model/experiment/multi_plane_geometry.h>
#include <dials/model/data/reflection.h>
#include "ray_predictor.h"

namespace dials { namespace algorithms {

  // Using lots of stuff from other namespaces
  using scitbx::vec2;
  using scitbx::vec3;
  using scitbx::mat3;
  using scitbx::af::double4;
  using dxtbx::model::Beam;
  using dxtbx::model::MultiPanelDetector;
  using dxtbx::model::Goniometer;
  using dxtbx::model::ScanData;
  using dials::model::BeamMultiPlaneIntersection;
  using dials::model::Reflection;

  // Typedef the miller_index type
  typedef cctbx::miller::index <> miller_index_type;

  /** A class to perform spot prediction. */
  class MultiPanelSpotPredictor {
  public:

    // A load of useful typedefs
    typedef Beam beam_type;
    typedef MultiPanelDetector detector_type;
    typedef Goniometer goniometer_type;
    typedef ScanData scan_type;
    typedef Reflection reflection_type;
    typedef MultiPanelDetector::coordinate_type coordinate_type;
    typedef scitbx::af::shared <reflection_type> reflection_list_type;

    /**
     * Initialise the spot predictor.
     * @param beam The beam parameters
     * @param detector The detector parameters
     * @param gonio The goniometer parameters
     * @param scan The scan parameters
     * @param unit_cell The unit cell parameters
     * @param space_group_type The space group struct
     * @param UB The ub matrix
     * @param d_min The resolution
     */
    MultiPanelSpotPredictor(const beam_type &beam,
                            const detector_type &detector,
                            const goniometer_type &gonio,
                            const scan_type &scan,
                            mat3 <double> UB)
      : predict_rays_(
          beam.get_direction(),
          gonio.get_rotation_axis(),
          UB,
          scan.get_oscillation_range()),
        plane_intersection_(
          detector.get_D_matrices(),
          get_image_extents(detector)),
        detector_(detector),
        scan_(scan) {}

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
    operator()(miller_index_type h) const {

      reflection_list_type reflections;

      // Predict the rays (phi, s1) for the miller index
      reflection_list_type rays = predict_rays_(h);

      // Loop through the reflections
      for (std::size_t i = 0; i < rays.size(); ++i) {

        // Get the ray data
        miller_index_type h = rays[i].get_miller_index();
        vec3<double> s1 = rays[i].get_beam_vector();
        double phi = rays[i].get_rotation_angle();

        // Try to calculate the detector coordinate in mm
        coordinate_type plane_xy_mm;
        try {
          plane_xy_mm = plane_intersection_(s1);
        } catch(error) {
          continue;
        }

        // Calculate the pixel coordinate of point
        int panel = plane_xy_mm.first;
        vec2<double> xy_mm = plane_xy_mm.second;
        vec2<double> xy_px = detector_.millimeter_to_pixel(plane_xy_mm);

        // Get the list of frames at which the reflection will be observed
        // and add the predicted observations to the list of reflections
        af::shared<double> frames = scan_.get_frames_with_angle(phi);
        for (std::size_t j = 0; j < frames.size(); ++j) {
          reflection_type r;
          r.set_miller_index(h);
          r.set_rotation_angle(phi);
          r.set_beam_vector(s1);
          r.set_image_coord_px(xy_px);
          r.set_image_coord_mm(xy_mm);
          r.set_frame_number(frames[j]);
          r.set_panel_number(panel);
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
    operator()(const const_ref<miller_index> &miller_indices) const {
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

    af::shared<double4> get_image_extents(const detector_type &detector) {
      af::shared<double4> extents(detector.num_panels(),
        af::init_functor_null<double4>());
      for (std::size_t i = 0; i < detector.num_panels(); ++i) {
        vec2<double> image_size_mm = detector[i].get_image_size_mm();
        extents[i] = double4(0.0, 0.0, image_size_mm[0], image_size_mm[1]);
      }
      return extents;
    }

    RayPredictor predict_rays_;
    BeamMultiPlaneIntersection plane_intersection_;
    const detector_type& detector_;
    const scan_type& scan_;
  };

}} // namespace dials::algorithms

#endif // DIALS_ALGORITHMS_SPOT_PREDICTION_MULTI_PANEL_SPOT_PREDICTOR_H

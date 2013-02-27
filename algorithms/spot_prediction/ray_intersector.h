/*
 * ray_intersector.h
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
#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/flex_types.h>
#include <cctbx/miller.h>
#include <dials/model/experiment/scan_helpers.h>
#include <dials/model/data/reflection.h>

namespace dials { namespace algorithms {

  // Using lots of stuff from other namespaces
  using scitbx::rad_as_deg;
  using scitbx::vec2;
  using scitbx::vec3;
  using scitbx::mat3;
  using model::mod_2pi;
  using model::is_angle_in_range;
  using model::Reflection;

  // Typedef the miller_index and flex_miller_index types
  typedef cctbx::miller::index <> miller_index;
  typedef scitbx::af::flex <miller_index> ::type flex_miller_index;

  /** A class to perform spot prediction. */
  template <typename DetectorType>
  class RayIntersector {
  public:

    // A load of useful typedefs
    typedef Reflection reflection_type;
    typedef scitbx::af::shared <reflection_type> reflection_list_type;

    /**
     * Initialise the ray predictor.
     * @param s0 The incident beam vector
     * @param m2 The rotation axis
     * @param UB The ub matrix
     * @param dphi The total oscillation range
     */
    RayIntersector() {}

    /** Virtual destructor to allow inheritance */
    virtual ~RayIntersector() {}

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
     * @param reflection The Reflection data
     * @returns The Reflection data with image volume coordinates
     */
    reflection_type
    operator()(const reflection_type &r) const {

      // Try to calculate the detector coordinate
      detector_coordinate_type coord;
      try {
        coord = get_detector_coord_(r.get_beam_vector());
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
      return r;
    }

    /**
     * For a given set of miller indices, predict the detector coordinates.
     * @param miller_indices The array of miller indices.
     */
    reflection_list_type
    operator()(const reflection_list_type &reflections) const {
      reflection_list_type reflections_new;
      for (std::size_t i = 0; i < reflections.size(); ++i) {
        try {
          reflection_type r = operator()(reflections[i]);
          reflections_new.push_back(r);
        } catch(error) {}
      }
      return reflections_new;
    }

  private:

    diffracted_beam_intersection <DetectorType> beam_coord_;
  };

  class RayIntersectorBase {};

  class RayIntersectorFlatPanel : public RayIntersectorBase {
  public:
    RayIntersectorFlatPanel() {}
    virtual ~RayIntersectorFlatPanel() {}

    virtual
    Reflection operator(const Reflection &r) {

    }

    virtual
    scitbx::af::shared <Reflection> operator(const scitbx::af::shared <Reflection> &r) {

    }

  private:
  };

  class RayIntersectorMultiFlatPanel : public RayIntersectorBase {
  public:
    RayIntersectorMultiFlatPanel() {}
    virtual ~RayIntersectorMultiFlatPanel() {}

    virtual
    ReflectionMultiPanel operator(const ReflectionMultiPanel &r) {

    }

    virtual
    scitbx::af::shared <ReflectionMultiPanel> operator(const scitbx::af::shared <ReflectionMultiPanel> &r) {

    }

  private:
  };

}} // namespace dials::algorithms

#endif // DIALS_ALGORITHMS_SPOT_PREDICTION_RAY_PREDICTOR_H

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
#ifndef DIALS_ALGORITHMS_SPOT_PREDICTION_RAY_INTERSECTOR_H
#define DIALS_ALGORITHMS_SPOT_PREDICTION_RAY_INTERSECTOR_H

#include <scitbx/vec2.h>
#include <scitbx/vec3.h>
#include <scitbx/array_family/flex_types.h>
#include <dxtbx/model/detector.h>
#include <dials/model/data/reflection.h>
#include <dials/error.h>

namespace dials { namespace algorithms {

  // Using lots of stuff from other namespaces
  using scitbx::vec2;
  using scitbx::vec3;
  using dxtbx::model::Detector;
  using dials::model::Reflection;
  using dials::model::ReflectionList;

  /** A class to perform spot prediction. */
  class RayIntersector {
  public:

    /**
     * Initialise the ray predictor.
     * @param Detector The detector
     */
    RayIntersector(const Detector &detector)
      : detector_(detector) {}

    /**
     * Calculate the spot locations on the detector image.
     * @param reflection The Reflection data
     * @returns The Reflection data with image volume coordinates
     */
    Reflection operator()(const Reflection &r) const {

      Reflection r_new(r);

      // Try to calculate the detector coordinate
      Detector::coord_type coord;
      try {
        coord = detector_.get_ray_intersection(r.get_beam_vector());
      } catch(dxtbx::error) {
        return r_new;
      }

      r_new.set_panel_number(coord.first);
      r_new.set_image_coord_mm(coord.second);

      return r_new;
    }

    Reflection operator()(const Reflection &r, int panel) const {

      Reflection r_new(r);

      // Try to calculate the detector coordinate
      vec2<double> coord;
      try {
        coord = detector_[panel].get_ray_intersection(r.get_beam_vector());
      } catch(dxtbx::error) {
        return r_new;
      }

      r_new.set_panel_number(panel);
      r_new.set_image_coord_mm(coord);

      return r_new;
    }


    /**
     * For a given set of miller indices, predict the detector coordinates.
     * @param miller_indices The array of miller indices.
     */
    ReflectionList
    operator()(const ReflectionList &reflections) const {
      ReflectionList reflections_new(reflections.size());
      for (std::size_t i = 0; i < reflections.size(); ++i) {
        Reflection r = operator()(reflections[i]);
        reflections_new[i] = r;
      }
      return reflections_new;
    }

    ReflectionList
    operator()(const ReflectionList &reflections, int panel) const {
      ReflectionList reflections_new(reflections.size());
      for (std::size_t i = 0; i < reflections.size(); ++i) {
        Reflection r = operator()(reflections[i], panel);
        reflections_new[i] = r;
      }
      return reflections_new;
    }

  private:

    Detector detector_;
  };

}} // namespace dials::algorithms

#endif // DIALS_ALGORITHMS_SPOT_PREDICTION_RAY_INTERSECTOR_H

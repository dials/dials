/*
 * ray_intersection.h
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
#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/flex_types.h>
#include <dxtbx/model/detector.h>
#include <dials/model/data/reflection.h>
#include <dials/error.h>

namespace dials { namespace algorithms {

  // Using lots of stuff from other namespaces
  using scitbx::vec2;
  using scitbx::vec3;
  using scitbx::af::shared;
  using dxtbx::model::Detector;
  using dials::model::Reflection;
  using dials::model::ReflectionList;

  /**
   * Calculate the intersection of a ray (given by the reflection object) with
   * the detector. Return as a new reflection object.
   * @param detector The detector model object
   * @param r The reflection object
   * @returns The new reflecton object
   */
  inline
  Reflection ray_intersection(const Detector &detector, const Reflection &r) {
    Reflection r_new(r);

    // Try to calculate the detector coordinate
    Detector::coord_type coord = detector.get_ray_intersection(
      r.get_beam_vector());

    // Set the panel and image coordinate number
    r_new.set_panel_number(coord.first);
    r_new.set_image_coord_mm(coord.second);
    r_new.set_image_coord_px(detector[coord.first].millimeter_to_pixel(
      coord.second));

    // Return reflection
    return r_new;
  }

  /**
   * Calculate the intersection of a ray (given by the reflection object) with
   * a specific detector panel. Return as a new reflection object.
   * @param detector The detector model object
   * @param r The reflection object
   * @param panel The panel number
   * @returns The new reflecton object
   */
  inline
  Reflection ray_intersection(const Detector &detector, const Reflection &r,
      std::size_t panel) {
    Reflection r_new(r);

    // Try to calculate the detector coordinate
    vec2<double> coord = detector[panel].get_ray_intersection(
      r.get_beam_vector());

    // Set the panel and image coordinate number
    r_new.set_panel_number(panel);
    r_new.set_image_coord_mm(coord);
    r_new.set_image_coord_px(detector[panel].millimeter_to_pixel(coord));

    // Return the reflection
    return r_new;
  }

  /**
   * Calculate the intersection of a list of ray (given by the reflection list)
   * with the detector. Return as a new reflection object.
   * @param detector The detector model object
   * @param r The reflection object
   * @returns The new reflecton object
   */
  inline
  shared<Reflection> ray_intersection(const Detector &detector,
      const ReflectionList &reflections) {
    shared<Reflection> reflections_new;
    reflections_new.reserve(reflections.size());
    for (std::size_t i = 0; i < reflections.size(); ++i) {
      try {
        reflections_new.push_back(ray_intersection(detector, reflections[i]));
      } catch(dxtbx::error) {
        // Do nothing
      } catch(dials::error) {
        // Do nothing
      }
    }
    return reflections_new;
  }

  /**
   * Calculate the intersection of a list of ray (given by the reflection list)
   * with a specific detector panel. Return as a new reflection object.
   * @param detector The detector model object
   * @param r The reflection object
   * @param panel The panel number
   * @returns The new reflecton object
   */
  inline
  shared<Reflection> ray_intersection(const Detector &detector,
      const ReflectionList &reflections, std::size_t panel) {
    shared<Reflection> reflections_new;
    reflections_new.reserve(reflections.size());
    for (std::size_t i = 0; i < reflections.size(); ++i) {
      try {
        reflections_new.push_back(ray_intersection(
          detector, reflections[i], panel));
      } catch(dxtbx::error) {
        // Do nothing
      }
    }
    return reflections_new;
  }

}} // namespace dials::algorithms

#endif // DIALS_ALGORITHMS_SPOT_PREDICTION_RAY_INTERSECTOR_H

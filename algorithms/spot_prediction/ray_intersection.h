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
#include <dxtbx/model/detector.h>
#include <dials/array_family/reflection_table.h>
#include <dials/model/data/reflection.h>
#include <dials/error.h>

namespace dials { namespace algorithms {

  // Using lots of stuff from other namespaces
  using scitbx::vec2;
  using scitbx::vec3;
  using dxtbx::model::Detector;
  using dxtbx::model::Panel;
  using dials::model::Reflection;


  //class ray_intersection2 {
  //public:

    //ray_intersection2(
        //const Detector &detector,
        //const af::const_ref< vec3<double> > &s1)
      //: xy_mm_(s1.size()),
        //xy_px_(s1.size()),
        //panel_(s1.size()) {

      //// Loop through and calculate the intersections
      //for (std::size_t i = 0; i < s1.size(); ++i) {
        //Detector::coord_type coord = detector.get_ray_intersection(s1[i]);
        //panel_[i] = coord.first;
        //xy_mm_[i] = coord.second;
        //xy_px_[i] = detector[panel_[i]].millimeter_to_pixel(xy_mm_[i]);
      //}
    //}

    //ray_intersection2(
        //const Detector &detector,
        //const af::const_ref< vec3<double> > &s1,
        //std::size_t panel)
      //: xy_mm_(s1.size()),
        //xy_px_(s1.size()),
        //panel_(s1.size(), panel) {

      //// Loop through and calculate the intersections
      //for (std::size_t i = 0; i < s1.size(); ++i) {
        //vec2<double> coord = detector[panel].get_ray_intersection(s1[i]);
        //xy_mm_[i] = coord;
        //xy_px_[i] = detector[panel_[i]].millimeter_to_pixel(xy_mm_[i]);
      //}
    //}

    //ray_intersection2(
        //const Detector &detector,
        //const af::const_ref< vec3<double> > &s1,
        //const af::const_ref< std::size_t > &panel)
      //: xy_mm_(s1.size()),
        //xy_px_(s1.size()),
        //panel_(panel) {

      //// Loop through and calculate the intersections
      //for (std::size_t i = 0; i < s1.size(); ++i) {
        //vec2<double> coord = detector[panel[i]].get_ray_intersection(s1[i]);
        //xy_mm_[i] = coord;
        //xy_px_[i] = detector[panel_[i]].millimeter_to_pixel(xy_mm_[i]);
      //}
    //}

    //af::shared< vec2<double> > xy_mm() const {
      //return xy_mm_;
    //}

    //af::shared< vec2<double> > xy_px() const {
      //return xy_px_;
    //}

    //af::shared< std::size_t > panel() const {
      //return panel_;
    //}

  //private:
    //af::shared< vec2<double> > xy_mm_;
    //af::shared< vec2<double> > xy_px_;
    //af::shared< std::size_t > panel_;
  //};


  struct Impact {
    vec3<double> px;
    vec3<double> mm;
  };

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
  af::shared<Reflection> ray_intersection(const Detector &detector,
      const af::const_ref<Reflection> &reflections) {
    af::shared<Reflection> reflections_new(af::reserve(reflections.size()));
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
  af::shared<Reflection> ray_intersection(const Detector &detector,
      const af::const_ref<Reflection> &reflections, std::size_t panel) {
    af::shared<Reflection> reflections_new(af::reserve(reflections.size()));
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

  inline
  af::shared<bool> ray_intersection(
      const Detector &detector,
      af::reflection_table reflections) {
    DIALS_ASSERT(reflections.is_consistent());
    DIALS_ASSERT(reflections.contains("s1"));
    DIALS_ASSERT(reflections.contains("phi"));
    af::const_ref< vec3<double> > s1 = reflections["s1"];
    af::const_ref< double > phi = reflections["phi"];
    af::ref<std::size_t> panel = reflections["panel"];
    af::ref< vec3<double> > xyzcalmm = reflections["xyzcal.mm"];
    af::shared<bool> success(reflections.size(), true);
    for (std::size_t i = 0; i < reflections.size(); ++i) {
      try {
      Detector::coord_type coord = detector.get_ray_intersection(s1[i]);
        xyzcalmm[i][0] = coord.second[0];
        xyzcalmm[i][1] = coord.second[1];
        xyzcalmm[i][2] = phi[i];
        panel[i] = coord.first;
      } catch(dxtbx::error) {
        success[i] = false;
      }
    }
    return success;
  }

  inline
  af::shared<bool> ray_intersection(
      const Detector &detector,
      af::reflection_table reflections,
      std::size_t panel_number) {
    DIALS_ASSERT(reflections.is_consistent());
    DIALS_ASSERT(reflections.contains("s1"));
    DIALS_ASSERT(reflections.contains("phi"));
    const Panel &p = detector[panel_number];
    af::const_ref< vec3<double> > s1 = reflections["s1"];
    af::const_ref< double > phi = reflections["phi"];
    af::ref<std::size_t> panel = reflections["panel"];
    af::ref< vec3<double> > xyzcalmm = reflections["xyzcal.mm"];
    af::shared<bool> success(reflections.size(), true);
    for (std::size_t i = 0; i < reflections.size(); ++i) {
      try {
        vec2<double> coord = p.get_ray_intersection(s1[i]);
        xyzcalmm[i][0] = coord[0];
        xyzcalmm[i][1] = coord[1];
        xyzcalmm[i][2] = phi[i];
        panel[i] = panel_number;
      } catch(dxtbx::error) {
        success[i] = false;
      }
    }
    return success;
  }

}} // namespace dials::algorithms

#endif // DIALS_ALGORITHMS_SPOT_PREDICTION_RAY_INTERSECTOR_H

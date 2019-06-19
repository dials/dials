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
#include <dials/error.h>

namespace dials { namespace algorithms {

  // Using lots of stuff from other namespaces
  using dxtbx::model::Detector;
  using dxtbx::model::Panel;
  using scitbx::vec2;
  using scitbx::vec3;

  // class ray_intersection2 {
  // public:

  // ray_intersection2(
  // const Detector &detector,
  // const af::const_ref< vec3<double> > &s1)
  //: xy_mm_(s1.size()),
  // xy_px_(s1.size()),
  // panel_(s1.size()) {

  //// Loop through and calculate the intersections
  // for (std::size_t i = 0; i < s1.size(); ++i) {
  // Detector::coord_type coord = detector.get_ray_intersection(s1[i]);
  // panel_[i] = coord.first;
  // xy_mm_[i] = coord.second;
  // xy_px_[i] = detector[panel_[i]].millimeter_to_pixel(xy_mm_[i]);
  //}
  //}

  // ray_intersection2(
  // const Detector &detector,
  // const af::const_ref< vec3<double> > &s1,
  // std::size_t panel)
  //: xy_mm_(s1.size()),
  // xy_px_(s1.size()),
  // panel_(s1.size(), panel) {

  //// Loop through and calculate the intersections
  // for (std::size_t i = 0; i < s1.size(); ++i) {
  // vec2<double> coord = detector[panel].get_ray_intersection(s1[i]);
  // xy_mm_[i] = coord;
  // xy_px_[i] = detector[panel_[i]].millimeter_to_pixel(xy_mm_[i]);
  //}
  //}

  // ray_intersection2(
  // const Detector &detector,
  // const af::const_ref< vec3<double> > &s1,
  // const af::const_ref< std::size_t > &panel)
  //: xy_mm_(s1.size()),
  // xy_px_(s1.size()),
  // panel_(panel) {

  //// Loop through and calculate the intersections
  // for (std::size_t i = 0; i < s1.size(); ++i) {
  // vec2<double> coord = detector[panel[i]].get_ray_intersection(s1[i]);
  // xy_mm_[i] = coord;
  // xy_px_[i] = detector[panel_[i]].millimeter_to_pixel(xy_mm_[i]);
  //}
  //}

  // af::shared< vec2<double> > xy_mm() const {
  // return xy_mm_;
  //}

  // af::shared< vec2<double> > xy_px() const {
  // return xy_px_;
  //}

  // af::shared< std::size_t > panel() const {
  // return panel_;
  //}

  // private:
  // af::shared< vec2<double> > xy_mm_;
  // af::shared< vec2<double> > xy_px_;
  // af::shared< std::size_t > panel_;
  //};

  inline af::shared<bool> ray_intersection(const Detector &detector,
                                           af::reflection_table reflections) {
    DIALS_ASSERT(reflections.is_consistent());
    DIALS_ASSERT(reflections.contains("s1"));
    DIALS_ASSERT(reflections.contains("phi"));
    af::const_ref<vec3<double> > s1 = reflections["s1"];
    af::const_ref<double> phi = reflections["phi"];
    af::ref<std::size_t> panel = reflections["panel"];
    af::ref<vec3<double> > xyzcalmm = reflections["xyzcal.mm"];
    af::shared<bool> success(reflections.size(), true);
    for (std::size_t i = 0; i < reflections.size(); ++i) {
      try {
        Detector::coord_type coord = detector.get_ray_intersection(s1[i]);
        xyzcalmm[i][0] = coord.second[0];
        xyzcalmm[i][1] = coord.second[1];
        xyzcalmm[i][2] = phi[i];
        panel[i] = coord.first;
      } catch (dxtbx::error) {
        success[i] = false;
      }
    }
    return success;
  }

  inline af::shared<bool> ray_intersection(const Detector &detector,
                                           af::reflection_table reflections,
                                           std::size_t panel_number) {
    DIALS_ASSERT(reflections.is_consistent());
    DIALS_ASSERT(reflections.contains("s1"));
    DIALS_ASSERT(reflections.contains("phi"));
    const Panel &p = detector[panel_number];
    af::const_ref<vec3<double> > s1 = reflections["s1"];
    af::const_ref<double> phi = reflections["phi"];
    af::ref<std::size_t> panel = reflections["panel"];
    af::ref<vec3<double> > xyzcalmm = reflections["xyzcal.mm"];
    af::shared<bool> success(reflections.size(), true);
    for (std::size_t i = 0; i < reflections.size(); ++i) {
      try {
        vec2<double> coord = p.get_ray_intersection(s1[i]);
        xyzcalmm[i][0] = coord[0];
        xyzcalmm[i][1] = coord[1];
        xyzcalmm[i][2] = phi[i];
        panel[i] = panel_number;
      } catch (dxtbx::error) {
        success[i] = false;
      }
    }
    return success;
  }

  inline af::shared<bool> ray_intersection(const Detector &detector,
                                           af::reflection_table reflections,
                                           af::ref<std::size_t> panel_numbers) {
    DIALS_ASSERT(reflections.is_consistent());
    DIALS_ASSERT(reflections.contains("s1"));
    DIALS_ASSERT(reflections.contains("phi"));
    DIALS_ASSERT(panel_numbers.size() == reflections.size());
    af::const_ref<vec3<double> > s1 = reflections["s1"];
    af::const_ref<double> phi = reflections["phi"];
    af::ref<std::size_t> panel = reflections["panel"];
    af::ref<vec3<double> > xyzcalmm = reflections["xyzcal.mm"];
    af::shared<bool> success(reflections.size(), true);
    for (std::size_t i = 0; i < reflections.size(); ++i) {
      try {
        const Panel &p = detector[panel_numbers[i]];
        vec2<double> coord = p.get_ray_intersection(s1[i]);
        xyzcalmm[i][0] = coord[0];
        xyzcalmm[i][1] = coord[1];
        xyzcalmm[i][2] = phi[i];
        panel[i] = panel_numbers[i];
      } catch (dxtbx::error) {
        success[i] = false;
      }
    }
    return success;
  }

}}  // namespace dials::algorithms

#endif  // DIALS_ALGORITHMS_SPOT_PREDICTION_RAY_INTERSECTOR_H

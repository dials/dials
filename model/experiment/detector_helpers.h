/*
 * detector_helpers.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_MODEL_EXPERIMENT_DETECTOR_HELPERS_H
#define DIALS_MODEL_EXPERIMENT_DETECTOR_HELPERS_H

//#include <boost/geometry.hpp>
//#include <boost/geometry/geometries/point.hpp>
//#include <boost/geometry/geometries/polygon.hpp>
#include <scitbx/vec2.h>
#include <scitbx/vec3.h>
#include <scitbx/array_family/tiny_types.h>

namespace dials { namespace model {

  using scitbx::vec2;
  using scitbx::vec3;
  using scitbx::af::double6;

  template <typename DetectorType>
  struct is_coordinate_valid {

    is_coordinate_valid(const DetectorType &detector)
      : image_size_(detector.get_image_size()) {}

    template <typename CoordinateType>
    bool operator()(CoordinateType coord) const {
      return (coord[0] >= 0 && coord[0] < image_size_[0])
          && (coord[1] >= 0 && coord[1] < image_size_[1]);
    }

  private:
    vec2 <std::size_t> image_size_;
  };

//  template <typename DetectorType>
//  struct mm_to_pixel {
//
//    mm_to_pixel(const DetectorType &detector) {}
//
//    vec2 <double> operator()(vec2 <double> mm) {
//
//    }
//  };

  
  template <typename DetectorType>
  struct diffracted_beam_detector_coord {

    diffracted_beam_detector_coord(const DetectorType &detector)
      : D_(detector.get_inverse_d_matrix()) {}

    vec2 <double> operator()(vec3 <double> s1) const {
      vec3 <double> v = D_ * s1;
      DIALS_ASSERT(v[2] > 0);
      return vec2 <double> (v[0] / v[2], v[1] / v[2]);
    }

  private:
    mat3 <double> D_;
  };

  /**
   * Get the image size in mm
   * @param detector The detector struct
   * @returns The detector image size in mm
   */
//  inline vec2 <double>
//  image_size_mm(const FlatPanelDetector &detector) {
//    return detector.get_image_size() * detector.get_pixel_size();
//  }

  /**
   * Get the pixel coordinate in mm in the laboratory frame
   * @param detector The detector struct
   * @param xy The xy pixel coordinate
   * @returns The detector pixel coordinate in mm in the laboratory frame
   */
//  template <typename T>
//  inline vec3 <double>
//  pixel_to_mm(const FlatPanelDetector &detector, vec2 <T> xy) {
//    return detector.get_d_matrix() * vec3 <double> (
//      (double) xy[0], (double) xy[1], 1.0);
//  }

  /**
   * Get the detector plane rectangle as lbx, lby, lbz, trx, try, trz
   * @param detector The detector struct
   * @returns The detector plane rectangle
   */
//  inline double6
//  plane_rectangle(const FlatPanelDetector &detector) {
//    vec3 <double> point1 = detector.get_origin();
//    vec3 <double> point2 = pixel_to_mm(detector, detector.get_image_size());
//    return double6(
//      point1[0], point1[1], point1[2],
//      point2[0], point2[1], point2[2]);
//  }

  /**
   * Check if the detector planes intersect.
   * @param a The first detector
   * @param b The second detector
   * @returns True/False do the detector planes intersect?
   */
//  inline bool
//  panels_intersect(const FlatPanelDetector &a, const FlatPanelDetector &b) {

//    using namespace boost::geometry;

//    typedef boost::geometry::model::point <double, 3, cs::cartesian> point;
//    typedef boost::geometry::model::polygon <point> polygon;

//    // Get the rectange of detector points
//    double6 rect_a = plane_rectangle(a);
//    double6 rect_b = plane_rectangle(b);

//    // Create a polygon for the panel a plane
//    polygon poly_a;
//    append(poly_a, point(rect_a[0], rect_a[1], rect_a[2]));
//    append(poly_a, point(rect_a[3], rect_a[1], rect_a[5]));
//    append(poly_a, point(rect_a[3], rect_a[4], rect_a[5]));
//    append(poly_a, point(rect_a[0], rect_a[4], rect_a[2]));
//    append(poly_a, point(rect_a[0], rect_a[1], rect_a[2]));

//    // Create a polygon for the panel b plane
//    polygon poly_b;
//    append(poly_b, point(rect_b[0], rect_b[1], rect_b[2]));
//    append(poly_b, point(rect_b[3], rect_b[1], rect_b[5]));
//    append(poly_b, point(rect_b[3], rect_b[4], rect_b[5]));
//    append(poly_b, point(rect_b[0], rect_b[4], rect_b[2]));
//    append(poly_b, point(rect_b[0], rect_b[1], rect_b[2]));

//    // Check if the polygons intersect
//    return intersects(poly_a, poly_b);
//  }

//  vec2 <double>
//  mm_to_pixel(const FlatPanelDetector &detector, vec3 <double> beam_vector) {
//    vec3 <double> v = detector.get_inverse_d_matrix() * beam_vector;
//    return vec2 <double> (v[0] / v[2], v[1] / v[2]);
//  }

  /**
   * Check if the detector coordinate is valid.
   * @param detector The MultiFlatPanelDetector struct
   * @returns True/False is the detector coordinate valid
   */
//  template <typename T>
//  inline bool
//  is_coordinate_valid(const MultiFlatPanelDetector &detector, vec3 <T> coord) {
//    int panel = (int)coord[0];
//    return (coord[0] >= 0 && coord[0] < detector.num_panels())
//        && (coord[1] >= 0 && coord[1] < detector[panel].get_image_size()[0])
//        && (coord[2] >= 0 && coord[2] < detector[panel].get_image_size()[1]);
//  }

}} // namespace dials::model

#endif // DIALS_MODEL_EXPERIMENT_DETECTOR_HELPERS_H

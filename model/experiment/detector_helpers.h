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
#include <scitbx/mat3.h>
#include <scitbx/array_family/tiny_types.h>
#include <scitbx/array_family/flex_types.h>
#include "detector.h"

namespace dials { namespace model {

  using scitbx::vec2;
  using scitbx::vec3;
  using scitbx::mat3;
  using scitbx::af::double6;

  // Create flex array typedefs
  typedef scitbx::af::flex<vec2<std::size_t> >::type flex_vec2_size_t;
  typedef scitbx::af::flex<mat3<double> >::type flex_mat3_double;

  /**
   * Declaration of a functor to check if a detector image coordinate is valid.
   * This functor is intended to be specialized for each detector type. For
   * instance a FlatPanelDetector object just needs to check an (x, y)
   * coordinate, whereas a MultiFlatPanelDetector object needs to check a
   * (panel, x, y) coordinate.
   * @tparam DetectorType The detector object type
   */
  template <typename DetectorType>
  struct is_coordinate_valid;

  /**
   * A functor to check if a FlatPanelDetector coordinate is valid. Check
   * that the (x, y) coordinate lies within the detector image size range.
   */
  template <>
  struct is_coordinate_valid <FlatPanelDetector> {

    /**
     * Initialise the functor with the image size.
     * @param detector The detector object
     */
    is_coordinate_valid(const FlatPanelDetector &detector)
      : image_size_(detector.get_image_size()) {}

    /**
     * Check the coordinate is valid.
     * @tparam CoordinateType The type of coordinate
     * @param coord The coordinate to check
     * @returns True/False the coordinate is valid
     */
    template <typename CoordinateType>
    bool operator()(CoordinateType coord) const {
      return (coord[0] >= 0 && coord[0] < image_size_[0])
          && (coord[1] >= 0 && coord[1] < image_size_[1]);
    }

  private:
    vec2 <std::size_t> image_size_;
  };

  /**
   * A functor to check if a MultiFlatPanelDetector coordinate is valid. Check
   * that the (panel x, y) coordinate lies within the number of panels and
   * the detector image size range.
   */
  template <>
  struct is_coordinate_valid <MultiFlatPanelDetector> {

    /**
     * Initialise the functor with the image size for each detector panel
     * @param detector The detector object
     */
    is_coordinate_valid(const MultiFlatPanelDetector &detector) 
      : image_size_(get_image_sizes(detector)),
        num_panels_(detector.num_panels()) {}

    /**
     * Check the coordinate is valid.
     * @tparam CoordinateType The type of coordinate
     * @param coord The coordinate to check
     * @returns True/False the coordinate is valid
     */
    template <typename CoordinateType>
    bool operator()(CoordinateType coord) const {
      int panel = coord[0]; 
      return (panel >= 0 && panel < num_panels_)
          && (coord[1] >= 0 && coord[1] < image_size_[panel][0])
          && (coord[2] >= 0 && coord[2] < image_size_[panel][1]);
    }

  private:
  
    /**
     * Put all the image sizes into a flex array
     * @param detector The detector object
     * @returns An array of image sizes.
     */
    flex_vec2_size_t get_image_sizes(
        const MultiFlatPanelDetector &detector) {
      flex_vec2_size_t result(detector.num_panels());
      for (std::size_t i = 0; i < result.size(); ++i) {
        result[i] = detector[i].get_image_size();
      }
      return result;
    }
  
    flex_vec2_size_t image_size_;
    std::size_t num_panels_;
  };

  /**
   * Declaration of a functor to to calculate the detector image coordinate
   * of a diffracted beam vector intersecting with the detector. This functor 
   * is intended to be specialized for each detector type.
   * @tparam DetectorType The type of the detector
   */
  template <typename DetectorType>
  struct diffracted_beam_to_pixel;
  
  /**
   * Specialization of the diffracted_beam_to_pixel functor for a generic
   * FlatPanelDetector object. Using a D matrix calculate the pixel coordinates
   */
  template <>
  struct diffracted_beam_to_pixel <FlatPanelDetector> {

    /** 
     * Initialise the transform from the detector object
     * @param detector The detector object
     */
    diffracted_beam_to_pixel(const FlatPanelDetector &detector)
      : D_(detector.get_inverse_d_matrix()),
        is_coord_valid_(detector) {}

    /**
     * Calculate the intersection point on the detector plane.
     * @param s1 The diffracted beam vector
     * @returns The (x, y) detector coordinate
     */
    vec2 <double> operator()(vec3 <double> s1) const {
      vec3 <double> v = D_ * s1;
      DIALS_ASSERT(v[2] > 0);
      vec2 <double> xy(v[0] / v[2], v[1] / v[2]);
      DIALS_ASSERT(is_coord_valid_(xy));
      return xy;
    }

  private:
    mat3 <double> D_;
    is_coordinate_valid <FlatPanelDetector> is_coord_valid_;
  };


  /**
   * Specialization of the diffracted_beam_to_pixel functor for a generic
   * MultiFlatPanelDetector object. Using a D matrix for each detextor panel,
   * calculate the pixel coordinates.
   * 
   * @todo This functor currently uses a brute force approach to finding the
   * detector panel which contains the reflection. All panels are checked
   * and where more than one panel could record the reflection, the closest
   * is chosen. This is ok for a small number of panels but doesn't scale
   * well for a large number.
   */
  template <>
  struct diffracted_beam_to_pixel <MultiFlatPanelDetector> {

    /** 
     * Initialise the transform from the detector object
     * @param detector The detector object
     */
    diffracted_beam_to_pixel(const MultiFlatPanelDetector &detector)
      : D_(get_inverse_d_matrices(detector)),
        is_coord_valid_(detector) {}

    /**
     * Calculate the intersection point on the detector planes.
     * @param s1 The diffracted beam vector
     * @returns The (panel, x, y) detector coordinate
     */
    vec3 <double> operator()(vec3 <double> s1) const {
      vec3 <double> pxy(-1, 0, 0);
      double w_max = 0;
      
      // Loop through all detectors. If the w component of the (u, v, w) 
      // vector points in the correct direction, then calculate the coordinate.
      // If the coordinate is valid and the w component is greater than that of
      // the current closest valid coordinate, then set this coordinate as the
      // current best bet. 
      for (std::size_t i = 0; i < D_.size(); ++i) {
        vec3 <double> v = D_[i] * s1;
        if (v[2] > 0) {
          vec3 <double> pxy_temp(i, v[0] / v[2], v[1] / v[2]);
          if (is_coord_valid_(pxy_temp) && v[2] > w_max) {
            pxy = pxy_temp;
            w_max = v[2];
          }
        }
      }
      
      // If no coordinate was found then raise an exception
      // otherwise return the coordinate.
      DIALS_ASSERT(w_max > 0);
      return pxy;
    }

  private:

    /**
     * Calculate the D matrix for each detector panel and cache.
     * @param detector The detector object
     * @returns An array of D matrices
     */
    flex_mat3_double get_inverse_d_matrices(
        const MultiFlatPanelDetector &detector) {
      flex_mat3_double result(detector.num_panels());
      for (std::size_t i = 0; i < result.size(); ++i) {
        result[i] = detector[i].get_inverse_d_matrix();
      }
      return result;
    }

    flex_mat3_double D_;
    is_coordinate_valid <MultiFlatPanelDetector> is_coord_valid_;
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

}} // namespace dials::model

#endif // DIALS_MODEL_EXPERIMENT_DETECTOR_HELPERS_H

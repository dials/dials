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

#include <scitbx/vec2.h>
#include <scitbx/vec3.h>
#include <scitbx/mat3.h>
#include <scitbx/array_family/tiny_types.h>
#include <scitbx/array_family/flex_types.h>
#include <dxtbx/model/detector.h>
#include <dxtbx/model/multi_panel_detector.h>
#include <dials/error.h>

namespace dials { namespace model {

  using scitbx::vec2;
  using scitbx::vec3;
  using scitbx::mat3;
  using scitbx::af::double6;
  using dxtbx::model::Detector;
  using dxtbx::model::MultiPanelDetector;

  // Create flex array typedefs
  typedef scitbx::af::flex<vec2<std::size_t> >::type flex_vec2_size_t;
  typedef scitbx::af::flex<vec2<double> >::type flex_vec2_double;
  typedef scitbx::af::flex<mat3<double> >::type flex_mat3_double;

  /**
   * Declaration of a functor to check if a detector image coordinate is valid.
   * This functor is intended to be specialized for each detector type. For
   * instance a Detector object just needs to check an (x, y)
   * coordinate, whereas a MultiPanelDetector object needs to check a
   * (panel, x, y) coordinate.
   * @tparam DetectorType The detector object type
   */
  template <typename DetectorType>
  struct is_coordinate_valid;

  /**
   * A functor to check if a Detector coordinate is valid. Check
   * that the (x, y) coordinate lies within the detector image size range.
   */
  template <>
  struct is_coordinate_valid <Detector> {

    /**
     * Initialise the functor with the image size.
     * @param detector The detector object
     */
    is_coordinate_valid(const Detector &detector)
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
   * A functor to check if a MultiPanelDetector coordinate is valid. Check
   * that the (panel x, y) coordinate lies within the number of panels and
   * the detector image size range.
   */
  template <>
  struct is_coordinate_valid <MultiPanelDetector> {

    /**
     * Initialise the functor with the image size for each detector panel
     * @param detector The detector object
     */
    is_coordinate_valid(const MultiPanelDetector &detector)
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
        const MultiPanelDetector &detector) {
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
  struct diffracted_beam_intersection;

  /**
   * Specialization of the diffracted_beam_to_pixel functor for a generic
   * Detector object. Using a D matrix calculate the pixel coordinates
   */
  template <>
  struct diffracted_beam_intersection <Detector> {

    /**
     * Initialise the transform from the detector object
     * @param detector The detector object
     */
    diffracted_beam_intersection(const Detector &detector)
      : D_(detector.get_D_matrix()),
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
    is_coordinate_valid <Detector> is_coord_valid_;
  };


  /**
   * Specialization of the diffracted_beam_to_pixel functor for a generic
   * MultiPanelDetector object. Using a D matrix for each detextor panel,
   * calculate the pixel coordinates.
   *
   * @todo This functor currently uses a brute force approach to finding the
   * detector panel which contains the reflection. All panels are checked
   * and where more than one panel could record the reflection, the closest
   * is chosen. This is ok for a small number of panels but doesn't scale
   * well for a large number.
   */
  template <>
  struct diffracted_beam_intersection <MultiPanelDetector> {

    /**
     * Initialise the transform from the detector object
     * @param detector The detector object
     */
    diffracted_beam_intersection(const MultiPanelDetector &detector)
      : D_(get_D_matrices(detector)),
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

    /**
     * Calculate the intersection point on a given detector plane.
     * @param s1 The diffracted beam vector
     * @param panel The panel
     * @returns The (panel, x, y) detector coordinate
     */
    vec3 <double> operator()(vec3 <double> s1, std::size_t panel) const {
      vec3 <double> v = D_[panel] * s1;
      DIALS_ASSERT(v[2] > 0);
      vec3 <double> pxy(panel, v[0] / v[2], v[1] / v[2]);
      return pxy;
    }

  private:

    /**
     * Calculate the D matrix for each detector panel and cache.
     * @param detector The detector object
     * @returns An array of D matrices
     */
    flex_mat3_double get_D_matrices(
        const MultiPanelDetector &detector) {
      flex_mat3_double result(detector.num_panels());
      for (std::size_t i = 0; i < result.size(); ++i) {
        result[i] = detector[i].get_D_matrix();
      }
      return result;
    }

    flex_mat3_double D_;
    is_coordinate_valid <MultiPanelDetector> is_coord_valid_;
  };

  /**
   * Template functor to map coordinates in the detector frame (in mm) to
   * pixel coordinates. This functor is specialized for each detector.
   * @tparam DetectorType The detector type.
   */
  template <typename DetectorType>
  struct millimeter_to_pixel;

  /**
   * Specialization to calculate millimerter to pixel mapping for generic
   * flat panel detector.
   */
  template <>
  struct millimeter_to_pixel <Detector> {

    /**
     * Initialise the functor with the image pixel size
     * @param detector The detector object
     */
    millimeter_to_pixel(const Detector &detector)
      : pixel_size_i_(1.0 / detector.get_pixel_size()) {}

    /**
     * Map the millimeter coordinate in the detector plane to pixels
     * @param xy The x, y detector coordinate in millimeters
     * @returns The x, y detector coordinate in pixels
     */
    vec2 <double> operator()(vec2 <double> xy) const {
      return xy * pixel_size_i_;
    }

  private:
    vec2 <double> pixel_size_i_;
  };

  /**
   * Specialization to calculate millimerter to pixel mapping for generic
   * multi flat panel detector.
   */
  template <>
  struct millimeter_to_pixel <MultiPanelDetector> {

    /**
     * Initialise the functor with the image pixel size
     * @param detector The detector object
     */
    millimeter_to_pixel(const MultiPanelDetector &detector)
      : pixel_size_i_(init_pixel_size_i(detector)) {}

    /**
     * Map the millimeter coordinate in the detector plane to pixels
     * @param pxy The panel, x, y detector coordinate in millimeters
     * @returns The panel, x, y detector coordinate in pixels
     */
    vec3 <double> operator()(vec3 <double> pxy) const {
      int panel = pxy[0];
      return vec3<double>(
        panel,
        pxy[1] * pixel_size_i_[panel][0],
        pxy[2] * pixel_size_i_[panel][1]);
    }

  private:

    /** Create an array of inverse pixel sizes */
    flex_vec2_double init_pixel_size_i(const MultiPanelDetector &detector) {
      flex_vec2_double result(detector.num_panels());
      for (std::size_t i = 0; i < detector.num_panels(); ++i) {
        result[i] = 1.0 / detector[i].get_pixel_size();
      }
      return result;
    }

    flex_vec2_double pixel_size_i_;
  };

  /**
   * Template functor to map coordinates in the detector frame (in pixels) to
   * millimeter coordinates. This functor is specialized for each detector.
   * @tparam DetectorType The detector type.
   */
  template <typename DetectorType>
  struct pixel_to_millimeter;

  /**
   * Specialization to calculate pixel to millimeter mapping for generic
   * flat panel detector.
   */
  template <>
  struct pixel_to_millimeter <Detector> {

    /**
     * Initialise the functor with the image pixel size
     * @param detector The detector object
     */
    pixel_to_millimeter(const Detector &detector)
      : pixel_size_(1.0 * detector.get_pixel_size()) {}

    /**
     * Map the pixel coordinate in the detector plane to millimeters
     * @param xy The x, y detector coordinate in pixels
     * @returns The x, y detector coordinate in millimeters
     */
    vec2 <double> operator()(vec2 <double> xy) const {
      return xy * pixel_size_;
    }

  private:
    vec2 <double> pixel_size_;
  };

  /**
   * Specialization to calculate pixel to millimeter mapping for generic
   * multi flat panel detector.
   */
  template <>
  struct pixel_to_millimeter <MultiPanelDetector> {

    /**
     * Initialise the functor with the image pixel size
     * @param detector The detector object
     */
    pixel_to_millimeter(const MultiPanelDetector &detector)
      : pixel_size_(init_pixel_size(detector)){}

    /**
     * Map the pixel coordinate in the detector plane to millimeters
     * @param pxy The panel, x, y detector coordinate in pixels
     * @returns The panel, x, y detector coordinate in millimeters
     */
    vec3 <double> operator()(vec3 <double> pxy) const {
      int panel = pxy[0];
      return vec3<double>(
        panel,
        pxy[1] * pixel_size_[panel][0],
        pxy[2] * pixel_size_[panel][1]);
    }

  private:

    /** Create an array of pixel sizes */
    flex_vec2_double init_pixel_size(const MultiPanelDetector &detector) {
      flex_vec2_double result(detector.num_panels());
      for (std::size_t i = 0; i < detector.num_panels(); ++i) {
        result[i] = 1.0 * detector[i].get_pixel_size();
      }
      return result;
    }

    flex_vec2_double pixel_size_;
  };

}} // namespace dials::model

#endif // DIALS_MODEL_EXPERIMENT_DETECTOR_HELPERS_H

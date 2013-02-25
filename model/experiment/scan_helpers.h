/*
 * scan_helpers.h
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#ifndef DIALS_MODEL_EXPERIMENT_SCAN_HELPERS_H
#define DIALS_MODEL_EXPERIMENT_SCAN_HELPERS_H

#include <cmath>
#include <scitbx/constants.h>
#include <scitbx/array_family/flex_types.h>
#include <dials/error.h>
#include "scan.h"

namespace dials { namespace model {

  using std::floor;
  using scitbx::vec2;
  using scitbx::af::flex_double;
  using scitbx::constants::two_pi;

  /** Convert the angle mod 2PI */
  inline
  double mod_2pi(double angle) {
    return angle - two_pi * floor(angle / two_pi);
  }

  /**
   * A functor to check if the angle is within the given range. The angular
   * range can be any two angles, plus or minus. The angle is to check can
   * also be any angle. The angle is considered within the range if the range
   * spans more than 2PI degrees and the angle is within the two range angles
   * when mod 2PI.
   */
  struct is_angle_in_range {

    /**
     * Initialise the functor, cache the range.
     * @param range The angular range
     */
    is_angle_in_range(vec2 <double> range)
      : range_(range) {}

    /**
     * Check the angle is within the range.
     * @param angle The angle to check
     * @returns True/False the angle is within the range
     */
    bool operator()(double angle) const {
      double diff_angle_range0 = mod_2pi(angle - range_[0]);
      double diff_angle_range1 = mod_2pi(angle - range_[1]);
      return range_[1] - range_[0] >= two_pi
          || diff_angle_range1 >= diff_angle_range0
          || diff_angle_range1 == 0;
    }

  private:
    vec2 <double> range_;
  };

  /**
   * A functor to get the range of equivalent angles mod 2PI degrees that lie
   * in the given angular range which can be more than 2PI degrees. For example
   * if the range is given as (A, B), the returned value will be (a, b) where
   * A <= a < A + 2PI; B - 2PI < b < B.
   */
  struct get_range_of_mod2pi_angles {

    /**
     * Initialise the functor with the given range. Range must be of the
     * form range[0] <= range[1].
     * @param range The angular range
     */
    get_range_of_mod2pi_angles(vec2 <double> range)
      : range_(range) {
      DIALS_ASSERT(range[0] <= range[1]);
    }

    /**
     * Calculate the range of angles.
     * @param angle The angle to check
     * @returns A pair of angles (a, b) where a = angle + n2PI and
     * b = angle + m2pi and both lie within the caches range values. If
     * b < a, the range is invalid.
     */
    vec2 <double> operator()(double angle) const {
      return vec2 <double> (
        angle - two_pi * floor((angle - range_[0]) / two_pi),
        angle + two_pi * floor((range_[1] - angle) / two_pi));
    }

  private:
    vec2 <double> range_;
  };

  /**
   * A functor to get the all the angles mod 2pi a given angle that lie
   * in the given angular range which can be more than 2pi degrees.
   */
  struct get_mod2pi_angles_in_range {

    /**
     * Initialise the functor with an angular range.
     * @param range The angular range.
     */
    get_mod2pi_angles_in_range(vec2 <double> range)
      : get_angular_range_(range) {}

    /**
     * Calculate and return an array of angles. If no angles are in the range,
     * then the array is empty.
     * @param angle The angle to use.
     * @returns An array of angles, a, where a = a + n2pi.
     */
    flex_double operator()(double angle) const {
      flex_double result;
      vec2 <double> range = get_angular_range_(angle);
      int n_angles = 1 + (int)floor((range[1] - range[0]) / two_pi);
      if (n_angles > 0) {
        result.resize(n_angles);
        for (std::size_t i = 0; i < n_angles; ++i) {
          result[i] = range[0] + i * two_pi;
        }
      }
      return result;
    }

  private:
    get_range_of_mod2pi_angles get_angular_range_;
  };

  /**
   * A functor to check if the angle is the range of angles coverd by the scan.
   * @tparam ScanType The type of scan object.
   */
  template <typename ScanType>
  struct is_scan_angle_valid;

  /**
   * A functor to check if the angle is the range of angles coverd by the scan.
   * The angle is checked using the is_angle_in_range object. The angular
   * range of the scan is given as (starting_angle, starting_angle + total
   * oscillation_range).
   * @tparam ScanType The type of scan object.
   */
  template <>
  struct is_scan_angle_valid <Scan> {

    /**
     * Initialise the functor. Use the scan data to initialise the base class.
     * @param scan The scan object.
     */
    is_scan_angle_valid(const Scan &scan)
      : is_angle_in_range_(vec2 <double> (
          scan.get_starting_angle(),
          scan.get_starting_angle() + scan.get_total_oscillation_range())) {}

    /**
     * Check the angle is within the range.
     * @param angle The angle to check
     * @returns True/False the angle is within the range
     */
    bool operator()(double angle) const {
      return is_angle_in_range_(angle);
    }

  private:
    is_angle_in_range is_angle_in_range_;
  };

  /**
   * A functor to calculate the angle corresponding to a scan frame.
   * @tparam ScanType The type of the scan object
   */
  template <typename ScanType>
  struct get_angle_from_frame;

  /**
   * A functor to calculate the angle corresponding to a scan frame.
   */
  template <>
  struct get_angle_from_frame <Scan> {

    /**
     * Initialise the functor with the scan object
     * @param scan The scan object.
     */
    get_angle_from_frame(const Scan &scan)
      : image_range_(scan.get_image_range()),
        starting_angle_(scan.get_starting_angle()),
        oscillation_range_(scan.get_oscillation_range()) {}

    /**
     * Calculate the angle corresponding to the given frame
     * @param frame The frame number
     * @returns The angle at the given frame
     */
    double operator()(double frame) const {
      return starting_angle_ + (frame - image_range_[0]) * oscillation_range_;
    }

  private:
    vec2 <int> image_range_;
    double starting_angle_;
    double oscillation_range_;
  };

  /**
   * A functor to calculate the frame corresponding to a scan angle. The raw
   * angle is used (i.e. not wrapped mod 2pi).
   * @tparam ScanType The type of the scan object
   */
  template <typename ScanType>
  struct get_frame_from_angle;

  /**
   * A functor to calculate the frame corresponding to a scan angle. The raw
   * angle is used (i.e. not wrapped mod 2pi).
   */
  template <>
  struct get_frame_from_angle <Scan> {

    /**
     * Initialise the functor with the scan object
     * @param scan The scan object.
     */
    get_frame_from_angle(const Scan &scan)
      : image_range_(scan.get_image_range()),
        starting_angle_(scan.get_starting_angle()),
        oscillation_range_(scan.get_oscillation_range()) {}

    /**
     * Calculate the frame corresponding to the given angle
     * @param angle The angle
     * @returns The frame at the given angle
     */
    double operator()(double angle) const {
      return image_range_[0] + (angle - starting_angle_) / oscillation_range_;
    }

  private:
    vec2 <int> image_range_;
    double starting_angle_;
    double oscillation_range_;
  };

  /**
   * A functor to calculate all the frames in the scan at which an observation
   * with a given angle will be observed. I.e. for a given angle, find all the
   * equivalent angles (i.e. mod 2pi) within the scan range and calculate the
   * frame number for each angle.
   * @tparam ScanType The type of the scan object
   */
  template <typename ScanType>
  struct get_all_frames_from_angle;

  /**
   * A functor to calculate all the frames in the scan at which an observation
   * with a given angle will be observed. I.e. for a given angle, find all the
   * equivalent angles (i.e. mod 2pi) within the scan range and calculate the
   * frame number for each angle.
   */
  template <>
  struct get_all_frames_from_angle <Scan> {

    /**
     * Initialise the functor with the scan object.
     * @param scan The scan object
     */
    get_all_frames_from_angle(const Scan &scan)
      : get_angles_(vec2 <double> (
          scan.get_starting_angle(),
          scan.get_starting_angle() + scan.get_total_oscillation_range())),
        get_frame_(scan) {}

    /**
     * Calculate and return an array of frame numbers at which a reflection
     * with a given rotation angle will be observed.
     * @param angle The rotation angle of the reflection
     * @returns The array of frame numbers
     */
    flex_double operator()(double angle) const {
      flex_double result = get_angles_(angle);
      for (std::size_t i = 0; i < result.size(); ++i) {
        result[i] = get_frame_(result[i]);
      }
      return result;
    }

  private:
    get_mod2pi_angles_in_range get_angles_;
    get_frame_from_angle <Scan> get_frame_;
  };

}} // namespace dials::model

#endif // DIALS_MODEL_EXPERIMENT_SCAN_HELPERS_H

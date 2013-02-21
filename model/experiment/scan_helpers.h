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
#include "scan.h"

namespace dials { namespace model {

  using std::floor;
  using scitbx::vec2;
  using scitbx::rad_as_deg;

  /** Convert the angle mod 360 */
  inline
  double mod_360(double angle) {
    return angle - 360.0 * floor(angle / 360);
  }

  /** Check if the angle is within the given range */
  template <bool deg = true>
  struct is_angle_in_range {

    /** Cache the range */
    is_angle_in_range(vec2 <double> range)
      : range_(deg ? range[0] : rad_as_deg(range[0]),
               deg ? range[1] : rad_as_deg(range[1])) {}

    /** Check the angle is within the range */
    bool operator()(double angle) const {
      if (!deg) {
        angle = rad_as_deg(angle);
      }
      double diff_angle_range0 = mod_360(angle - range_[0]);
      double diff_angle_range1 = mod_360(angle - range_[1]);
      return range_[1] - range_[0] >= 360.0
          || diff_angle_range1 >= diff_angle_range0
          || diff_angle_range1 == 0;
    }

  private:
    vec2 <double> range_;
  };

  /** Check the angle is valid */
  template <typename ScanType, bool deg = true>
  struct is_scan_angle_valid
    : is_angle_in_range <deg> {

    /** Set the angular range from the starting angle and oscillation range */
    is_scan_angle_valid(const ScanType &scan)
      : is_angle_in_range <deg> (vec2 <double> (
          scan.get_starting_angle(),
          scan.get_starting_angle() + scan.get_total_oscillation_range())) {}
  };

}} // namespace dials::model

#endif // DIALS_MODEL_EXPERIMENT_SCAN_HELPERS_H

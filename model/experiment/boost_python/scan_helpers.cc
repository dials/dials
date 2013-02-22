/*
 * scan_helpers.cc
 *
 *  Copyright (C) 2013 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <boost/format.hpp>
#include <string>
#include <scitbx/constants.h>
#include <dials/model/experiment/scan.h>
#include <dials/model/experiment/scan_helpers.h>

namespace dials { namespace model { namespace boost_python {

  using namespace boost::python;

  using scitbx::deg_as_rad;

  vec2 <double> deg_as_rad(vec2 <double> angles) {
    return vec2 <double> (deg_as_rad(angles[0]), deg_as_rad(angles[1]));
  }

  bool is_angle_in_range_wrapper(vec2 <double> range, double angle, bool deg) {
    return is_angle_in_range(deg ? deg_as_rad(range) : range)(
        deg ? deg_as_rad(angle) : angle);
  }

  bool is_scan_angle_valid_wrapper(const Scan &scan, double angle, bool deg) {
    return is_scan_angle_valid<Scan>(scan)(deg ? deg_as_rad(angle) : angle);
  }

  void export_scan_helpers()
  {
    def("is_angle_in_range", 
      &is_angle_in_range_wrapper, (
        arg("range"),
        arg("angle"),
        arg("deg") = false));
    def("is_scan_angle_valid", 
      &is_scan_angle_valid_wrapper, (
        arg("scan"),
        arg("angle"),
        arg("deg") = false));
  }

}}} // namespace = dials::model::boost_python

/*
 * detector_helpers.cc
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
#include <dials/model/experiment/detector.h>
#include <dials/model/experiment/detector_helpers.h>

namespace dials { namespace model { namespace boost_python {

  using namespace boost::python;

  template <typename DetectorType, typename CoordinateType>
  bool is_coordinate_valid_wrapper(const DetectorType &detector, 
      CoordinateType coord) {
    return is_coordinate_valid <DetectorType>(detector)(coord);
  }

  void export_detector_helpers()
  {
    def("is_coordinate_valid", 
      &is_coordinate_valid_wrapper <FlatPanelDetector, vec2 <double> >, (
        arg("detector"), 
        arg("coordinate")));
    def("is_coordinate_valid", 
      &is_coordinate_valid_wrapper <FlatPanelDetector, vec2 <int> >, (
        arg("detector"),
        arg("coordinate")));
  }

}}} // namespace = dials::model::boost_python

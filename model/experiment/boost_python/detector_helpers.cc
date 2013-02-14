/*
 * scan.cc
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
#include <dials/model/experiment/detector_helpers.h>

using namespace boost::python;

namespace dials { namespace model { namespace experiment { namespace detector { namespace boost_python {

  void export_detector_helpers()
  {
    bool (*flat_panel_detector_is_coordinate_valid)(
      const FlatPanelDetector &, vec2 <double>) = &is_coordinate_valid;

    bool (*multi_flat_panel_detector_is_coordinate_valid)(
      const MultiFlatPanelDetector &, vec3 <double>) = &is_coordinate_valid;

    def("is_coordinate_valid", flat_panel_detector_is_coordinate_valid);
    def("is_coordinate_valid", mulit_flat_panel_detector_is_coordinate_valid);
  }

}}}} // namespace = dials::model::experiment::boost_python

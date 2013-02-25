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

  void export_flat_panel_detector_helpers()
  {
//    // Export the is_coordinate_valid functor
//    class_<is_coordinate_valid <FlatPanelDetector> >("is_coord_valid", no_init)
//      .def(init<const FlatPanelDetector&>())
//      .def("__call__", &is_coordinate_valid <FlatPanelDetector>::operator()<vec2<int> >)
//      .def("__call__", &is_coordinate_valid <FlatPanelDetector>::operator()<vec2<double> >);
  }

  void export_multi_flat_panel_detector_helpers()
  {
    // Export the is_coordinate_valid functor  
//    class_<is_coordinate_valid <MultiFlatPanelDetector> >("is_coord_valid.MultiFlatPanelDetector")
//      .def(init<const MultiFlatPanelDetector&>())
//      .def("__call__", &MultiFlatPanelDetector::operator()<vec3<int> >)
//      .def("__call__", &MultiFlatPanelDetector::operator()<vec3<double> >);
  }

  void export_detector_helpers()
  {
    export_flat_panel_detector_helpers();
    export_multi_flat_panel_detector_helpers();
  }

}}} // namespace = dials::model::boost_python

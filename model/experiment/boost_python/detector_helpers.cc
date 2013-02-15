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

namespace dials { namespace model { namespace boost_python {

  using namespace boost::python;

  void export_detector_helpers()
  {
    bool (*flat_panel_detector_is_coordinate_valid)(
      const FlatPanelDetector &, vec2 <double>) = &is_coordinate_valid <double>;

    bool (*multi_flat_panel_detector_is_coordinate_valid)(
      const MultiFlatPanelDetector &, vec3 <double>) = &is_coordinate_valid <double>;

    def("is_coordinate_valid", flat_panel_detector_is_coordinate_valid);
    def("is_coordinate_valid", multi_flat_panel_detector_is_coordinate_valid);
    def("image_size_mm", &image_size_mm);
    def("pixel_to_mm", &pixel_to_mm <double>);
    def("plane_rectangle", &plane_rectangle);
    def("panels_intersect", &panels_intersect);
  }

}}} // namespace = dials::model::boost_python

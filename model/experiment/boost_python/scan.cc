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
#include <dials/model/experiment/scan.h>

namespace dials { namespace model { namespace boost_python {

  using namespace boost::python;

  std::string scan_to_string(const Scan &scan) {
    boost::format fmt(
      "Scan:\n"
      "    image range:       (%1%, %2%)\n"
      "    starting angle:    %3%\n"
      "    oscillation range: %4%");
        
    fmt % scan.get_image_range()[0];
    fmt % scan.get_image_range()[1];
    fmt % scan.get_starting_angle();
    fmt % scan.get_oscillation_range();
    return fmt.str();
  }

  void export_scan()
  {
    // Export ScanBase
    class_ <ScanBase> ("ScanBase");

    // Export Scan : ScanBase
    class_ <Scan, bases <ScanBase> > ("Scan")
      .def(init <vec2 <int>, double, double> ((
          arg("image_range"), 
          arg("starting_angle"),
          arg("oscillation_range"))))
      .add_property("image_range",  
        &Scan::get_image_range,
        &Scan::set_image_range)
      .add_property("starting_angle",  
        &Scan::get_starting_angle,
        &Scan::set_starting_angle)
      .add_property("oscillation_range",  
        &Scan::get_oscillation_range,
        &Scan::set_oscillation_range)
      .add_property("total_oscillation_range",
        &Scan::get_total_oscillation_range)
      .add_property("num_images",
        &Scan::get_num_images)
      .def("__eq__", &Scan::operator==)
      .def("__nq__", &Scan::operator!=)
      .def("__str__", &scan_to_string);
  }

}}} // namespace = dials::model::boost_python

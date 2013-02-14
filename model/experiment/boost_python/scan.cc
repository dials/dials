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

using namespace boost::python;

namespace dials { namespace model { namespace experiment { namespace boost_python {

  std::string scan_to_string(const Scan &scan) {
    boost::format fmt(
      "Scan:\n"
      "    starting frame:    %1%\n"
      "    starting angle:    %2%\n"
      "    oscillation range: %3%\n"
      "    num frames:        %4%");
        
    fmt % scan.get_starting_frame();
    fmt % scan.get_starting_angle();
    fmt % scan.get_oscillation_range();
    fmt % scan.get_num_frames();
    return fmt.str();
  }

  void export_scan()
  {
    // Export ScanBase
    class_ <ScanBase> ("ScanBase");

    // Export Scan : ScanBase
    class_ <Scan, bases <ScanBase> > ("Scan")
      .def(init <int, double, double, int> ((
          arg("starting_frame"), 
          arg("starting_angle"),
          arg("oscillation_range"),
          arg("num_frames"))))
      .add_property("starting_frame",  
        &Scan::get_starting_frame,
        &Scan::set_starting_frame)
      .add_property("starting_angle",  
        &Scan::get_starting_angle,
        &Scan::set_starting_angle)
      .add_property("oscillation_range",  
        &Scan::get_oscillation_range,
        &Scan::set_oscillation_range)
      .add_property("num_frames",  
        &Scan::get_num_frames,
        &Scan::set_num_frames)
      .add_property("total_oscillation_range",
        &Scan::get_total_oscillation_range)
      .def("__str__", &scan_to_string);
  }

}}}} // namespace = dials::model::experiment::boost_python

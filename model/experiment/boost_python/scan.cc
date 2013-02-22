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
#include <boost/python/make_constructor.hpp>
#include <boost/format.hpp>
#include <string>
#include <scitbx/constants.h>
#include <dials/model/experiment/scan.h>

namespace dials { namespace model { namespace boost_python {

  using namespace boost::python;
  using scitbx::deg_as_rad;
  using scitbx::rad_as_deg;

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

  static Scan* make_scan(vec2 <int> image_range, double starting_angle,
      double oscillation_range, bool deg) {
    Scan *scan = NULL;
    if (deg) {
      scan = new Scan(image_range, deg_as_rad(starting_angle), 
        deg_as_rad(oscillation_range));
    } else {
      scan = new Scan(image_range, starting_angle, oscillation_range);
    }
    return scan;
  }
  
  static double get_starting_angle_deg(Scan &scan) {
    return rad_as_deg(scan.get_starting_angle());
  }

  static void set_starting_angle_deg(Scan &scan, double angle) {
    scan.set_starting_angle(deg_as_rad(angle));
  }

  static double get_oscillation_range_deg(Scan &scan) {
    return rad_as_deg(scan.get_oscillation_range());
  }

  static void set_oscillation_range_deg(Scan &scan, double angle) {
    scan.set_oscillation_range(deg_as_rad(angle));
  }

  static double get_total_oscillation_range_deg(Scan &scan) {
    return rad_as_deg(scan.get_total_oscillation_range());
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
      .def("__init__", 
        make_constructor(
          &make_scan, 
          default_call_policies(), (
          arg("image_range"),
          arg("starting_angle"),
          arg("oscillation_range"),
          arg("deg"))))
      .add_property("image_range",  
        &Scan::get_image_range,
        &Scan::set_image_range)
      .add_property("starting_angle",  
        &Scan::get_starting_angle,
        &Scan::set_starting_angle)
      .add_property("starting_angle_deg",
        &get_starting_angle_deg,
        &set_starting_angle_deg)
      .add_property("oscillation_range",  
        &Scan::get_oscillation_range,
        &Scan::set_oscillation_range)
      .add_property("oscillation_range_deg",  
        &get_oscillation_range_deg,
        &set_oscillation_range_deg)
      .add_property("total_oscillation_range",
        &Scan::get_total_oscillation_range)
      .add_property("total_oscillation_range_deg",
        &get_total_oscillation_range_deg)
      .add_property("num_images",
        &Scan::get_num_images)
      .def("__eq__", &Scan::operator==)
      .def("__nq__", &Scan::operator!=)
      .def("__str__", &scan_to_string);
  }

}}} // namespace = dials::model::boost_python

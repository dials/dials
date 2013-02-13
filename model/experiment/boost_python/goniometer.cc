/*
 * goniometer.cc
 *
 *   Copyright (C) 2013 Diamond Light Source, James Parkhurst
 *
 *   This code is distributed under the BSD license, a copy of which is
 *   included in the root directory of this package.
 */
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <boost/format.hpp>
#include <string>
#include <dials/model/experiment/goniometer.h>

using namespace boost::python;

namespace dials { namespace model { namespace experiment { namespace boost_python {

  std::string goniometer_to_string(const Goniometer &goniometer) {
    boost::format fmt(
      "Goniometer:\n"
      "    rotation axis:  (%1%, %2%, %3%)\n"
      "    fixed rotation: [%4%, %5%, %6%,\n"
      "                     %7%, %8%, %9%,\n"
      "                     %10%, %11%, %12%]");        
    fmt % goniometer.get_rotation_axis()[0];
    fmt % goniometer.get_rotation_axis()[1];
    fmt % goniometer.get_rotation_axis()[2];
    fmt % goniometer.get_fixed_rotation_matrix()[0];
    fmt % goniometer.get_fixed_rotation_matrix()[1];
    fmt % goniometer.get_fixed_rotation_matrix()[2];
    fmt % goniometer.get_fixed_rotation_matrix()[3];
    fmt % goniometer.get_fixed_rotation_matrix()[4];
    fmt % goniometer.get_fixed_rotation_matrix()[5];
    fmt % goniometer.get_fixed_rotation_matrix()[6];
    fmt % goniometer.get_fixed_rotation_matrix()[7];
    fmt % goniometer.get_fixed_rotation_matrix()[8];
    return fmt.str();
  }

  std::string kappa_goniometer_direction_to_string(
      KappaGoniometer::Direction d) {
    if (d == KappaGoniometer::PlusY) {
      return "+y";
    } else if (d == KappaGoniometer::PlusZ) {
      return "+z";
    } else if (d == KappaGoniometer::MinusY) {
      return "-y";
    } else if (d == KappaGoniometer::MinusZ) {
      return "-z";
    }
    return "none";
  }

  std::string kappa_goniometer_scan_axis_to_string(
      KappaGoniometer::ScanAxis a) {
    if (a == KappaGoniometer::Omega) {
      return "omega";
    } else if (a == KappaGoniometer::Phi) {
      return "phi";
    }
    return "none";
  }

  std::string kappa_goniometer_to_string(const KappaGoniometer &goniometer) {
    boost::format fmt(
      "%1%\n"
      "    alpha angle:    %2%\n"
      "    omega angle:    %3%\n"
      "    kappa angle:    %4%\n"
      "    phi angle:      %5%\n"
      "    direction:      %6%\n"
      "    scan axis:      %7%\n"
      "    omega axis:     (%8%, %9%, %10%)\n"
      "    kappa axis:     (%11%, %12%, %13%)\n"
      "    phi axis:       (%14%, %15%, %16%)");
    fmt % goniometer_to_string(goniometer);
    fmt % goniometer.get_alpha_angle();
    fmt % goniometer.get_omega_angle();
    fmt % goniometer.get_kappa_angle();
    fmt % goniometer.get_phi_angle();
    fmt % kappa_goniometer_direction_to_string(goniometer.get_direction());
    fmt % kappa_goniometer_scan_axis_to_string(goniometer.get_scan_axis());
    fmt % goniometer.get_omega_axis()[0];
    fmt % goniometer.get_omega_axis()[1];
    fmt % goniometer.get_omega_axis()[2];
    fmt % goniometer.get_kappa_axis()[0];
    fmt % goniometer.get_kappa_axis()[1];
    fmt % goniometer.get_kappa_axis()[2];
    fmt % goniometer.get_phi_axis()[0];
    fmt % goniometer.get_phi_axis()[1];
    fmt % goniometer.get_phi_axis()[2];
    return fmt.str();
  }

  void export_goniometer() 
  {
    class_ <GoniometerBase> ("GoniometerBase");

    class_ <Goniometer, bases <GoniometerBase> > ("Goniometer")
      .def(init <vec3 <double> > ((
          arg("rotation_axis"))))
      .def(init <vec3 <double>,
                 mat3 <double> > ((
          arg("rotation_axis"), 
          arg("fixed_rotation_matrix"))))
      .add_property("rotation_axis",  
        &Goniometer::get_rotation_axis,
        &Goniometer::set_rotation_axis)
      .add_property("fixed_rotation_matrix",  
        &Goniometer::get_fixed_rotation_matrix,
        &Goniometer::set_fixed_rotation_matrix)
      .def("__str__", &goniometer_to_string);
     
    enum_ <KappaGoniometer::Direction> ("KappaDirection")
      .value("PlusY", KappaGoniometer::PlusY)
      .value("PlusZ", KappaGoniometer::PlusZ)
      .value("MinusY", KappaGoniometer::MinusY)
      .value("MinusZ", KappaGoniometer::MinusZ);

    enum_ <KappaGoniometer::ScanAxis> ("KappaScanAxis")
      .value("Omega", KappaGoniometer::Omega)
      .value("Phi", KappaGoniometer::Phi);

    class_ <KappaGoniometer, bases <Goniometer> > ("KappaGoniometer")
      .def(init <double, 
                 double, 
                 double, 
                 double, 
                 KappaGoniometer::Direction, 
                 KappaGoniometer::ScanAxis> ((
          arg("alpha"),
          arg("omega"),
          arg("kappa"),
          arg("phi"),
          arg("direction"),
          arg("scan_axis"))))
      .add_property("alpha_angle",
        &KappaGoniometer::get_alpha_angle)
      .add_property("omega_angle",
        &KappaGoniometer::get_omega_angle)
      .add_property("kappa_angle",
        &KappaGoniometer::get_kappa_angle)
      .add_property("phi_angle",
        &KappaGoniometer::get_phi_angle)
      .add_property("direction",
        &KappaGoniometer::get_direction)
      .add_property("scan_axis",
        &KappaGoniometer::get_scan_axis)
      .add_property("omega_axis",
        &KappaGoniometer::get_omega_axis)
      .add_property("kappa_axis",
        &KappaGoniometer::get_kappa_axis)
      .add_property("phi_axis",
        &KappaGoniometer::get_phi_axis)
      .def("__str__", &kappa_goniometer_to_string);
  }

}}}} // namespace dials::model::experiment::boost_python

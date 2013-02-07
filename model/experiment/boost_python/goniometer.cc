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

  void export_goniometer() 
  {
    class_ <Goniometer> ("Goniometer")
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
  }

}}}} // namespace dials::model::experiment::boost_python

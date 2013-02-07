/*
 * beam.cc
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
#include <dials/model/experiment/beam.h>

using namespace boost::python;

namespace dials { namespace model { namespace experiment { namespace boost_python {

  std::string beam_to_string(const Beam &beam) {
    boost::format fmt(
      "Beam:\n"
      "    direction :            (%1%, %2%, %3%)\n"
      "    wavelength:            %4%\n"
      "    polarization:          (%5%, %6%, %7%)\n"
      "    polarization fraction: %8%");
        
    fmt % beam.get_direction()[0];
    fmt % beam.get_direction()[1];
    fmt % beam.get_direction()[2];
    fmt % beam.get_wavelength();
    fmt % beam.get_polarization()[0];
    fmt % beam.get_polarization()[1];
    fmt % beam.get_polarization()[2];
    fmt % beam.get_polarization_fraction();
    return fmt.str();
  }

  void export_beam()
  {
    class_ <Beam> ("Beam")
      .def(init <double> ((
          arg("wavelength"))))
      .def(init <double,
                 vec3 <double> > ((
          arg("wavelength"), 
          arg("direction"))))
      .def(init <double,
                 vec3 <double>, 
                 vec3 <double>, 
                 double> ((
          arg("wavelength"),
          arg("direction"), 
          arg("polarization"),
          arg("polarization_fraction"))))
      .def(init <vec3 <double>, 
                 vec3 <double>, 
                 double> ((
          arg("direction"), 
          arg("polarization"),
          arg("polarization_fraction"))))
      .def(init <vec3 <double> > ((
          arg("direction"))))
      .def(init <vec3 <double>, 
                 vec3 <double>, 
                 double> ((
          arg("direction"), 
          arg("polarization"),
          arg("polarization_fraction"))))
      .add_property("direction",  
        &Beam::get_direction,
        &Beam::set_direction)
      .add_property("wavelength", 
        &Beam::get_wavelength)
      .add_property("polarization",
        &Beam::get_polarization,
        &Beam::set_polarization)
      .add_property("polarization_fraction",
        &Beam::get_polarization_fraction,
        &Beam::set_polarization_fraction)
      .def("__str__", &beam_to_string);
  }

}}}} // namespace = dials::model::experiment::boost_python

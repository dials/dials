/*
 * reflection.cc
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
#include <scitbx/array_family/flex_types.h>
#include <scitbx/array_family/boost_python/flex_wrapper.h>
#include <dials/model/data/reflection.h>

namespace dials { namespace model { namespace boost_python {

  using namespace boost::python;

  std::string reflection_base_to_string(const ReflectionBase &reflection) {
    boost::format fmt(
      "Reflection:\n"
      "    miller index:    (%1%, %2%, %3%)");
        
    fmt % reflection.get_miller_index()[0];
    fmt % reflection.get_miller_index()[1];
    fmt % reflection.get_miller_index()[2];
    return fmt.str();
  }

  std::string reflection_to_string(
      const Reflection &reflection) {
    boost::format fmt(
      "Reflection:\n"
      "    miller index:     (%1%, %2%, %3%)\n"
      "    rotation angle:   %4%\n"
      "    beam vector:      (%5%, %6%, %7%)\n"
      "    image coord (px): (%8%, %9%)\n"
      "    image coord (mm): (%10%, %11%)\n"
      "    frame number:     %12%\n"
      "    panel number:     %13%");
        
    fmt % reflection.get_miller_index()[0];
    fmt % reflection.get_miller_index()[1];
    fmt % reflection.get_miller_index()[2];
    fmt % reflection.get_rotation_angle();
    fmt % reflection.get_beam_vector()[0];
    fmt % reflection.get_beam_vector()[1];
    fmt % reflection.get_beam_vector()[2];
    fmt % reflection.get_image_coord_px()[0];
    fmt % reflection.get_image_coord_px()[1];
    fmt % reflection.get_image_coord_mm()[0];
    fmt % reflection.get_image_coord_mm()[1];
    fmt % reflection.get_frame_number();
    fmt % reflection.get_panel_number();
    return fmt.str();
  }

  void export_reflection()
  {
    class_<ReflectionBase>("ReflectionBase")
      .def(init <miller_index_type> ((
          arg("miller_index"))))
      .add_property("miller_index", 
        &Reflection::get_miller_index,
        &Reflection::set_miller_index)
      .def("is_zero", &Reflection::is_zero)
      .def("__str__", &reflection_base_to_string);

    class_<Reflection, bases<ReflectionBase> > ("Reflection")
      .def(init <const Reflection &>())
      .def(init <miller_index_type> ((
          arg("miller_index"))))
      .add_property("rotation_angle", 
        &Reflection::get_rotation_angle,
        &Reflection::set_rotation_angle)
      .add_property("beam_vector", 
        &Reflection::get_beam_vector,
        &Reflection::set_beam_vector)
      .add_property("image_coord_mm",
        &Reflection::get_image_coord_mm,
        &Reflection::set_image_coord_mm)
      .add_property("image_coord_px",
        &Reflection::get_image_coord_px,
        &Reflection::set_image_coord_px)
      .add_property("frame_number",
        &Reflection::get_frame_number,
        &Reflection::set_frame_number)
      .add_property("panel_number",
        &Reflection::get_panel_number,
        &Reflection::set_panel_number)      
      .add_property("shoebox",
        &Reflection::get_shoebox,
        &Reflection::set_shoebox)  
      .def("__str__", &reflection_to_string);          

    scitbx::af::boost_python::flex_wrapper <Reflection>::plain("ReflectionList");        
  }

}}} // namespace dials::model::boost_python

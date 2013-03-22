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
      "    panel number:     %13%\n"
      "    shoebox:          (%14%, %15%, %16%, %17%, %18%, %19%)\n");
//      "    centroid pos:     (%20%, %21%, %22%)\n"
//      "    centroid var:     (%23%, %24%, %24%)\n");
        
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
    fmt % reflection.get_shoebox()[0];
    fmt % reflection.get_shoebox()[1];
    fmt % reflection.get_shoebox()[2];
    fmt % reflection.get_shoebox()[3];
    fmt % reflection.get_shoebox()[4];
    fmt % reflection.get_shoebox()[5];
//    fmt % reflection.get_centroid_position()[0];
//    fmt % reflection.get_centroid_position()[1];
//    fmt % reflection.get_centroid_position()[2];
//    fmt % reflection.get_centroid_variance()[0];
//    fmt % reflection.get_centroid_variance()[1];
//    fmt % reflection.get_centroid_variance()[2];
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
      .add_property("image",
        &Reflection::get_image,
        &Reflection::set_image)
      .add_property("image_weights",
        &Reflection::get_image_weights,
        &Reflection::set_image_weights)
      .add_property("transformed_image",
        &Reflection::get_transformed_image,
        &Reflection::set_transformed_image)
      .add_property("centroid_position",
        &Reflection::get_centroid_position,
        &Reflection::set_centroid_position)
      .add_property("centroid_variance",
        &Reflection::get_centroid_variance,
        &Reflection::set_centroid_variance)
      .def("__str__", &reflection_to_string);          

    scitbx::af::boost_python::flex_wrapper <Reflection>::plain("ReflectionList");        
  }

}}} // namespace dials::model::boost_python

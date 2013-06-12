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
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <sstream>
#include <string>
#include <scitbx/array_family/flex_types.h>
#include <scitbx/array_family/boost_python/flex_wrapper.h>
#include <scitbx/array_family/boost_python/flex_pickle_double_buffered.h>
#include <dials/model/data/reflection.h>

namespace dials { namespace model { namespace boost_python {

  using namespace boost::python;

  using scitbx::af::boost_python::flex_pickle_double_buffered;

  std::string reflection_to_string(const Reflection &reflection) {
    std::stringstream ss;
    ss << reflection;
    return ss.str();
  }

  struct ReflectionPickleSuite : boost::python::pickle_suite {
//    static
//    boost::python::tuple getinitargs(const Reflection &r) {
//      // ...
//    }

    static
    boost::python::tuple getstate(boost::python::object obj) {
      const Reflection &r = extract<const Reflection&>(obj)();
      return boost::python::make_tuple(
        obj.attr("__dict__"),
        r.get_miller_index(),
        r.get_rotation_angle(),
        r.get_beam_vector(),
        r.get_image_coord_mm(),
        r.get_image_coord_px(),
        r.get_frame_number(),
        r.get_panel_number(),
        r.get_bounding_box(),
        r.get_centroid_position(),
        r.get_centroid_variance(),
        r.get_centroid_sq_width(),
        r.get_shoebox(),
        r.get_shoebox_mask(),
        r.get_transformed_shoebox());
    }

    static
    void setstate(boost::python::object obj, boost::python::tuple state) {
      Reflection &r = extract<Reflection&>(obj)();
      
      // Check that the number of items is correct
      if (len(state) != 15) {
        PyErr_SetObject(PyExc_ValueError, (
          "expected 15-item tuple in call to __setstate__; got %s" 
          % state).ptr());
        throw_error_already_set();
      }

      // restore the object's __dict__
      dict d = extract<dict>(obj.attr("__dict__"))();
      d.update(state[0]);
      
      // restore the internal state of the C++ reflection object
      r.set_miller_index(extract<cctbx::miller::index<> >(state[1]));
      r.set_rotation_angle(extract<double>(state[2]));
      r.set_beam_vector(extract<vec3<double> >(state[3]));
      r.set_image_coord_mm(extract<vec2<double> >(state[4]));
      r.set_image_coord_px(extract<vec2<double> >(state[5]));
      r.set_frame_number(extract<double>(state[6]));
      r.set_panel_number(extract<int>(state[7]));
      r.set_bounding_box(extract<int6>(state[8]));
      r.set_centroid_position(extract<vec3<double> >(state[9]));
      r.set_centroid_variance(extract<vec3<double> >(state[10]));
      r.set_centroid_sq_width(extract<vec3<double> >(state[11]));
      r.set_shoebox(extract<const flex_int&>(state[12]));
      r.set_shoebox_mask(extract<const flex_int&>(state[13]));
      r.set_transformed_shoebox(extract<const flex_double&>(state[14]));
    }
  };

  void export_reflection()
  {
    class_<ReflectionBase>("ReflectionBase")
      .def(init <miller_index_type> ((
          arg("miller_index"))))
      .add_property("miller_index", 
        &ReflectionBase::get_miller_index,
        &ReflectionBase::set_miller_index)
      .add_property("entering",
        &ReflectionBase::get_entering,
        &ReflectionBase::set_entering)
      .add_property("status",
        &ReflectionBase::get_status,
        &ReflectionBase::set_status)
      .def("is_valid", &ReflectionBase::is_valid)
      .def("is_active", &ReflectionBase::is_active)
      .def("set_valid", &ReflectionBase::set_valid)
      .def("set_active", &ReflectionBase::set_active)
      .def("is_zero", &ReflectionBase::is_zero);

    flex_int& (Reflection::*reflection_get_shoebox)() = 
      &Reflection::get_shoebox;

    flex_int& (Reflection::*reflection_get_shoebox_mask)() = 
      &Reflection::get_shoebox_mask;

    flex_int& (Reflection::*reflection_get_shoebox_background)() = 
      &Reflection::get_shoebox_background;

    flex_double& (Reflection::*reflection_get_transformed_shoebox)() = 
      &Reflection::get_transformed_shoebox;

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
      .add_property("bounding_box",
        &Reflection::get_bounding_box,
        &Reflection::set_bounding_box)  
      .add_property("shoebox",
        make_function(reflection_get_shoebox, 
          return_internal_reference<>()),
        &Reflection::set_shoebox)
      .add_property("shoebox_mask",
        make_function(reflection_get_shoebox_mask, 
          return_internal_reference<>()),
        &Reflection::set_shoebox_mask)
      .add_property("shoebox_background",
        make_function(reflection_get_shoebox_background, 
          return_internal_reference<>()),
        &Reflection::set_shoebox_background)        
      .add_property("transformed_shoebox",
        make_function(reflection_get_transformed_shoebox, 
          return_internal_reference<>()),
        &Reflection::set_transformed_shoebox)
      .add_property("centroid_position",
        &Reflection::get_centroid_position,
        &Reflection::set_centroid_position)
      .add_property("centroid_variance",
        &Reflection::get_centroid_variance,
        &Reflection::set_centroid_variance)
      .add_property("centroid_sq_width",
        &Reflection::get_centroid_sq_width,
        &Reflection::set_centroid_sq_width)
      .add_property("intensity",
        &Reflection::get_intensity,
        &Reflection::set_intensity)
      .add_property("intensity_variance",
        &Reflection::get_intensity_variance,
        &Reflection::set_intensity_variance)
      .def("__str__", &reflection_to_string)
      .def_pickle(ReflectionPickleSuite());          

//      class_<ReflectionList>("ReflectionList")
//        .def(vector_indexing_suite<ReflectionList>())
//        .enable_pickling();
    scitbx::af::boost_python::flex_wrapper 
      <Reflection, return_internal_reference<> >::plain("ReflectionList")
        .enable_pickling();        
  }

}}} // namespace dials::model::boost_python

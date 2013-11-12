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
#include <scitbx/array_family/ref_reductions.h>
#include <dials/model/data/reflection.h>
#include <dials/error.h>
#include <dials/model/data/boost_python/pickle.h>
#include <dials/model/data/observation.h>
#include <dials/model/data/shoebox.h>
#include <scitbx/array_family/boost_python/flex_pickle_double_buffered.h>

namespace dials { namespace model { namespace boost_python {

  using namespace boost::python;
  using scitbx::af::flex_double;
  using scitbx::af::flex_int;
  using scitbx::af::flex;

  typedef cctbx::miller::index<> miller_index;

  std::string reflection_to_string(const Reflection &reflection) {
    std::stringstream ss;
    ss << reflection;
    return ss.str();
  }

  af::flex<Reflection>::type* init_from_observation_and_shoebox(
      const af::const_ref<Observation> &o, 
      const af::const_ref< Shoebox<Reflection::float_type> > &s) {
    DIALS_ASSERT(o.size() == s.size());
    af::shared<Reflection> result(o.size());   
    for (std::size_t i = 0; i < result.size(); ++i) {

      // Check panel numbers
      DIALS_ASSERT(o[i].panel == s[i].panel);
      result[i].panel_number_ = o[i].panel;

      // Copy observation info
      result[i].centroid_position_ = o[i].centroid.px.position;
      result[i].centroid_variance_ = o[i].centroid.px.std_err_sq;
      result[i].centroid_sq_width_ = o[i].centroid.px.variance;
      result[i].intensity_ = o[i].intensity.observed.value;
      result[i].intensity_variance_ = o[i].intensity.observed.variance;
      result[i].corrected_intensity_ = o[i].intensity.corrected.value;
      result[i].corrected_intensity_variance_ = o[i].intensity.corrected.variance;
      
      // Copy shoebox info
      result[i].bounding_box_ = s[i].bbox;
      result[i].shoebox_ = s[i].data;
      result[i].shoebox_mask_ = s[i].mask;
      result[i].shoebox_background_ = s[i].background;
    }
    
    return new af::flex<Reflection>::type(
      result, af::flex_grid<>(result.size()));
  }
  
  
  void set_shoebox(Reflection &obj, 
      flex<Reflection::float_type>::type &data) {
    DIALS_ASSERT(data.accessor().all().size() == 3);
    obj.set_shoebox(af::versa<Reflection::float_type, af::c_grid<3> >(
      data.handle(), af::c_grid<3>(data.accessor())));
  }

//  static
//  flex_double get_shoebox(Reflection &obj) {
//    return flex_double(obj.get_shoebox().handle(), 
//      obj.get_shoebox().accessor().as_flex_grid());
//  }
  
  void set_shoebox_mask(Reflection &obj, flex_int &data) {
    DIALS_ASSERT(data.accessor().all().size() == 3);
    obj.set_shoebox_mask(af::versa<int, af::c_grid<3> >(
      data.handle(), af::c_grid<3>(data.accessor())));
  }

//  static
//  flex_int get_shoebox_mask(Reflection &obj) {
//    return flex_int(obj.get_shoebox_mask().handle(), 
//      obj.get_shoebox_mask().accessor().as_flex_grid());
//  }
  
  void set_shoebox_background(Reflection &obj, 
      flex<Reflection::float_type>::type &data) {
    DIALS_ASSERT(data.accessor().all().size() == 3);
    obj.set_shoebox_background(af::versa<Reflection::float_type, af::c_grid<3> >(
      data.handle(), af::c_grid<3>(data.accessor())));
  }

//  static
//  flex_double get_shoebox_background(Reflection &obj) {
//    return flex_double(obj.get_shoebox_background().handle(), 
//      obj.get_shoebox_background().accessor().as_flex_grid());
//  }
  
  void set_transformed_shoebox(Reflection &obj, 
      flex<Reflection::float_type>::type &data) {
    DIALS_ASSERT(data.accessor().all().size() == 3);
    obj.set_transformed_shoebox(af::versa<Reflection::float_type, af::c_grid<3> >(
      data.handle(), af::c_grid<3>(data.accessor())));
  }

//  static
//  flex_double get_transformed_shoebox(Reflection &obj) {
//    return flex_double(obj.get_transformed_shoebox().handle(), 
//      obj.get_transformed_shoebox().accessor().as_flex_grid());
//  }
  
  void set_transformed_shoebox_background(Reflection &obj, 
      flex<Reflection::float_type>::type &data) {
    DIALS_ASSERT(data.accessor().all().size() == 3);
    obj.set_transformed_shoebox_background(af::versa<Reflection::float_type, af::c_grid<3> >(
      data.handle(), af::c_grid<3>(data.accessor())));
  }

  static
  af::shared< cctbx::miller::index<> > get_miller_index(const af::const_ref<Reflection> &r) {
    af::shared< cctbx::miller::index<> > result(r.size());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = r[i].get_miller_index();
    }
    return result;
  }

  static
  af::shared< double > get_rotation_angle(const af::const_ref<Reflection> &r) {
    af::shared< double > result(r.size());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = r[i].get_rotation_angle();
    }
    return result;
  }

  static
  af::shared< vec3<double> > get_beam_vector(const af::const_ref<Reflection> &r) {
    af::shared< vec3<double> > result(r.size());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = r[i].get_beam_vector();
    }
    return result;
  }

  static
  void set_beam_vector(
      af::ref<Reflection> const& r,
      af::const_ref<vec3 <double> > const& beam_vector) {
    DIALS_ASSERT(r.size() == beam_vector.size());
    for (std::size_t i = 0; i < r.size(); ++i) {
      r[i].set_beam_vector(beam_vector[i]);
    }
  }

  static
  af::shared< vec3<double> > get_centroid_position(const af::const_ref<Reflection> &r) {
    af::shared< vec3<double> > result(r.size());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = r[i].get_centroid_position();
    }
    return result;
  }

  static
  af::shared< vec2<double> > get_image_coord_px(const af::const_ref<Reflection> &r) {
    af::shared< vec2<double> > result(r.size());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = r[i].get_image_coord_px();
    }
    return result;
  }
  
  static
  af::shared< vec2<double> > get_image_coord_mm(const af::const_ref<Reflection> &r) {
    af::shared< vec2<double> > result(r.size());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = r[i].get_image_coord_mm();
    }
    return result;
  }
  
  static
  af::shared< double > get_frame_number(const af::const_ref<Reflection> &r) {
    af::shared< double > result(r.size());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = r[i].get_frame_number();
    }
    return result;
  }

  static
  af::shared< int > get_crystal(const af::const_ref<Reflection> &r) {
    af::shared< int > result(r.size());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = r[i].get_crystal();
    }
    return result;
  }

  static
  void set_crystal(const af::ref<Reflection> &r,
                                   const af::const_ref<int> &crystal) {
    DIALS_ASSERT(r.size() == crystal.size());
    for (std::size_t i = 0; i < r.size(); ++i) {
      r[i].set_crystal(crystal[i]);
    }
  }

  static
  af::shared< bool > get_is_valid(const af::const_ref<Reflection> &r) {
    af::shared< bool > result(r.size());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = r[i].is_valid();
    }
    return result;
  } 
  
//  static
//  flex_double get_transformed_shoebox_background(Reflection &obj) {
//    return flex_double(obj.get_transformed_shoebox_background().handle(), 
//      obj.get_transformed_shoebox_background().accessor().as_flex_grid());
//  }
  
  void reflection_wrapper()
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
      .def("is_strong", &ReflectionBase::is_strong)
      .def("set_valid", &ReflectionBase::set_valid)
      .def("set_active", &ReflectionBase::set_active)
      .def("set_strong", &ReflectionBase::set_strong)
      .def("is_zero", &ReflectionBase::is_zero);

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
//      .def("shoebox", &get_shoebox)
//      .def("shoebox", &set_shoebox)
//      .def("shoebox_mask", &get_shoebox_mask)
//      .def("shoebox_mask", &set_shoebox_mask)
//      .def("shoebox_background", &get_shoebox_background)
//      .def("shoebox_background", &set_shoebox_background)      
//      .def("transformed_shoebox", &get_transformed_shoebox)
//      .def("transformed_shoebox", &set_transformed_shoebox)   
//      .def("transformed_shoebox_background", &get_transformed_shoebox_background)
//      .def("transformed_shoebox_background", &set_transformed_shoebox_background)      
      .add_property("shoebox",
        &Reflection::get_shoebox,
        &set_shoebox)
      .add_property("shoebox_mask",
        &Reflection::get_shoebox_mask,
        &set_shoebox_mask)
      .add_property("shoebox_background",
        &Reflection::get_shoebox_background,
        &set_shoebox_background)        
      .add_property("transformed_shoebox",
        &Reflection::get_transformed_shoebox,
        &set_transformed_shoebox)
      .add_property("transformed_shoebox_background",
        &Reflection::get_transformed_shoebox_background,
        &set_transformed_shoebox_background)        
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
      .add_property("corrected_intensity",
        &Reflection::get_corrected_intensity,
        &Reflection::set_corrected_intensity)
      .add_property("corrected_intensity_variance",
        &Reflection::get_corrected_intensity_variance,
        &Reflection::set_corrected_intensity_variance)        
      .add_property("crystal",
        &Reflection::get_crystal,
        &Reflection::set_crystal)        
      .def("__str__", &reflection_to_string)
      .def_pickle(reflection::ReflectionPickleSuite());          

    scitbx::af::boost_python::flex_wrapper 
      <Reflection, return_internal_reference<> >::plain("ReflectionList")
        .def("__init__", make_constructor(
          &init_from_observation_and_shoebox))
        .def_pickle(scitbx::af::boost_python::flex_pickle_double_buffered<
          Reflection, 
          reflection::to_string, 
          reflection::from_string>())
        .def("miller_index", &get_miller_index)
        .def("rotation_angle", &get_rotation_angle)
        .def("centroid_position", &get_centroid_position)
        .def("beam_vector", &get_beam_vector)
        .def("set_beam_vector", &set_beam_vector)
        .def("crystal", &get_crystal)
        .def("set_crystal", &set_crystal)
        .def("image_coord_px", &get_image_coord_px)
        .def("image_coord_mm", &get_image_coord_mm)
        .def("frame_number", &get_frame_number)
        .def("is_valid", &get_is_valid);
  }
  
  void export_reflection() {
    reflection_wrapper();
  }

}}} // namespace dials::model::boost_python

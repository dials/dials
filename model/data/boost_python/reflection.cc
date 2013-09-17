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

  typedef cctbx::miller::index<> miller_index;

  std::string reflection_to_string(const Reflection &reflection) {
    std::stringstream ss;
    ss << reflection;
    return ss.str();
  }

  af::flex<Reflection>::type* init_from_observation_and_shoebox(
      const af::const_ref<Observation> &o, 
      const af::const_ref<Shoebox> &s) {
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
    
    return new af::flex<Reflection>::type(result, af::flex_grid<>(result.size()));
  }
  
  af::shared<miller_index> reflection_list_get_miller_index(const af::const_ref<Reflection> &r) {
    af::shared<miller_index> result(r.size());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = r[i].get_miller_index();
    }
    return result;
  }
  
  af::shared<int> reflection_list_get_status(const af::const_ref<Reflection> &r) {
    af::shared<int> result(r.size());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = r[i].get_status();
    }
    return result;
  }
  
  af::shared<bool> reflection_list_get_entering(const af::const_ref<Reflection> &r) {
    af::shared<bool> result(r.size());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = r[i].get_entering();
    }
    return result;
  }

  af::shared<double> reflection_list_get_rotation_angle(const af::const_ref<Reflection> &r) {
    af::shared<double> result(r.size());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = r[i].get_rotation_angle();
    }
    return result;
  }

  af::shared< vec3<double> > reflection_list_get_beam_vector(const af::const_ref<Reflection> &r) {
    af::shared< vec3<double> > result(r.size());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = r[i].get_beam_vector();
    }
    return result;
  }
  
  af::shared< vec2<double> > reflection_list_get_image_coord_mm(const af::const_ref<Reflection> &r) {
    af::shared< vec2<double> > result(r.size());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = r[i].get_image_coord_mm();
    }
    return result;
  }
  
  af::shared< vec2<double> > reflection_list_get_image_coord_px(const af::const_ref<Reflection> &r) {
    af::shared< vec2<double> > result(r.size());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = r[i].get_image_coord_px();
    }
    return result;
  }
  
  af::shared<double> reflection_list_get_frame_number(const af::const_ref<Reflection> &r) {
    af::shared<double> result(r.size());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = r[i].get_frame_number();
    }
    return result;
  }
  
  af::shared<int> reflection_list_get_panel_number(const af::const_ref<Reflection> &r) {
    af::shared<int> result(r.size());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = r[i].get_panel_number();
    }
    return result;
  }
  
  af::shared<int6> reflection_list_get_bounding_box(const af::const_ref<Reflection> &r) {
    af::shared<int6> result(r.size());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = r[i].get_bounding_box();
    }
    return result;
  }
  
  af::shared< vec3<double> > reflection_list_get_centroid_position(const af::const_ref<Reflection> &r) {
    af::shared< vec3<double> > result(r.size());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = r[i].get_centroid_position();
    }
    return result;
  }
  
  af::shared< vec3<double> > reflection_list_get_centroid_variance(const af::const_ref<Reflection> &r) {
    af::shared< vec3<double> > result(r.size());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = r[i].get_centroid_variance();
    }
    return result;
  }
  
  af::shared< vec3<double> > reflection_list_get_centroid_sq_width(const af::const_ref<Reflection> &r) {
    af::shared< vec3<double> > result(r.size());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = r[i].get_centroid_sq_width();
    }
    return result;
  }
  
  af::shared<double> reflection_list_get_intensity(const af::const_ref<Reflection> &r) {
    af::shared<double> result(r.size());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = r[i].get_intensity();
    }
    return result;
  }
  
  af::shared<double> reflection_list_get_intensity_variance(const af::const_ref<Reflection> &r) {
    af::shared<double> result(r.size());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = r[i].get_intensity_variance();
    }
    return result;
  }  
  
  af::shared<double> reflection_list_get_corrected_intensity(const af::const_ref<Reflection> &r) {
    af::shared<double> result(r.size());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = r[i].get_corrected_intensity();
    }
    return result;
  }  
 
  af::shared<double> reflection_list_get_corrected_intensity_variance(
      const af::const_ref<Reflection> &r) {
    af::shared<double> result(r.size());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = r[i].get_corrected_intensity_variance();
    }
    return result;
  }   
 
  af::shared<double> reflection_list_get_shoebox(const af::const_ref<Reflection> &r) {
    std::size_t result_size = 0;
    for (std::size_t i = 0; i < r.size(); ++i) {
      result_size += r[i].get_shoebox().size();
    }
    af::shared<double> result(result_size);
    for (std::size_t i = 0, k = 0; i < r.size(); ++i) {
      af::versa< double, af::c_grid<3> > s = r[i].get_shoebox();
      for (std::size_t j = 0; j < s.size(); ++j) {
        result[k++] = s[j];
      }
    }
    return result;
  }
  
  af::shared<int> reflection_list_get_shoebox_mask(const af::const_ref<Reflection> &r) {
    std::size_t result_size = 0;
    for (std::size_t i = 0; i < r.size(); ++i) {
      result_size += r[i].get_shoebox_mask().size();
    }
    af::shared<int> result(result_size);
    for (std::size_t i = 0, k = 0; i < r.size(); ++i) {
      af::versa<int, af::c_grid<3> > s = r[i].get_shoebox_mask();
      for (std::size_t j = 0; j < s.size(); ++j) {
        result[k++] = s[j];
      }
    }
    return result;
  }
  
  af::shared<double> reflection_list_get_shoebox_background(const af::const_ref<Reflection> &r) {
    std::size_t result_size = 0;
    for (std::size_t i = 0; i < r.size(); ++i) {
      result_size += r[i].get_shoebox_background().size();
    }
    af::shared<double> result(result_size);
    for (std::size_t i = 0, k = 0; i < r.size(); ++i) {
      af::versa< double, af::c_grid<3> > s = r[i].get_shoebox_background();
      for (std::size_t j = 0; j < s.size(); ++j) {
        result[k++] = s[j];
      }
    }
    return result;
  }
  
  af::shared<double> reflection_list_get_transformed_shoebox(const af::const_ref<Reflection> &r) {
    std::size_t result_size = 0;
    for (std::size_t i = 0; i < r.size(); ++i) {
      result_size += r[i].get_transformed_shoebox().size();
    }
    af::shared<double> result(result_size);
    for (std::size_t i = 0, k = 0; i < r.size(); ++i) {
      af::versa< double, af::c_grid<3> > s = r[i].get_transformed_shoebox();
      for (std::size_t j = 0; j < s.size(); ++j) {
        result[k++] = s[j];
      }
    }
    return result;
  }
  
  
  void reflection_list_set_miller_index(af::ref<Reflection> r, 
      const af::const_ref<miller_index> &input) {
    DIALS_ASSERT(input.size() == r.size());        
    for (std::size_t i = 0; i < input.size(); ++i) {
      r[i].set_miller_index(input[i]);
    }
  }

  void reflection_list_set_status(af::ref<Reflection> r, 
      const af::const_ref<int> &input) {
    DIALS_ASSERT(input.size() == r.size());        
    for (std::size_t i = 0; i < input.size(); ++i) {
      r[i].set_status(input[i]);
    }
  }

  void reflection_list_set_entering(af::ref<Reflection> r, 
      const af::const_ref<bool> &input) {
    DIALS_ASSERT(input.size() == r.size());        
    for (std::size_t i = 0; i < input.size(); ++i) {
      r[i].set_entering(input[i]);
    }
  }


  void reflection_list_set_rotation_angle(af::ref<Reflection> r, 
      const af::const_ref<double> &input) {
    DIALS_ASSERT(input.size() == r.size());        
    for (std::size_t i = 0; i < input.size(); ++i) {
      r[i].set_rotation_angle(input[i]);
    }
  }

  void reflection_list_set_beam_vector(af::ref<Reflection> r, 
      const af::const_ref< vec3<double> > &input) {
    DIALS_ASSERT(input.size() == r.size());        
    for (std::size_t i = 0; i < input.size(); ++i) {
      r[i].set_beam_vector(input[i]);
    }
  }

  void reflection_list_set_image_coord_mm(af::ref<Reflection> r, 
      const af::const_ref< vec2<double> > &input) {
    DIALS_ASSERT(input.size() == r.size());        
    for (std::size_t i = 0; i < input.size(); ++i) {
      r[i].set_image_coord_mm(input[i]);
    }
  }

  void reflection_list_set_image_coord_px(af::ref<Reflection> r, 
      const af::const_ref< vec2<double> > &input) {
    DIALS_ASSERT(input.size() == r.size());        
    for (std::size_t i = 0; i < input.size(); ++i) {
      r[i].set_image_coord_px(input[i]);
    }
  }
  
  void reflection_list_set_frame_number(af::ref<Reflection> r, 
      const af::const_ref<double> &input) {
    DIALS_ASSERT(input.size() == r.size());        
    for (std::size_t i = 0; i < input.size(); ++i) {
      r[i].set_frame_number(input[i]);
    }
  }

  void reflection_list_set_panel_number(af::ref<Reflection> r, 
      const af::const_ref<int> &input) {
    DIALS_ASSERT(input.size() == r.size());        
    for (std::size_t i = 0; i < input.size(); ++i) {
      r[i].set_panel_number(input[i]);
    }
  }

  void reflection_list_set_bounding_box(af::ref<Reflection> r, 
      const af::const_ref<int6> &input) {
    DIALS_ASSERT(input.size() == r.size());        
    for (std::size_t i = 0; i < input.size(); ++i) {
      r[i].set_bounding_box(input[i]);
    }
  }
  
  void reflection_list_set_centroid_position(af::ref<Reflection> r, 
      const af::const_ref< vec3<double> > &input) {
    DIALS_ASSERT(input.size() == r.size());        
    for (std::size_t i = 0; i < input.size(); ++i) {
      r[i].set_centroid_position(input[i]);
    }
  }
  
  void reflection_list_set_centroid_variance(af::ref<Reflection> r, 
      const af::const_ref< vec3<double> > &input) {
    DIALS_ASSERT(input.size() == r.size());        
    for (std::size_t i = 0; i < input.size(); ++i) {
      r[i].set_centroid_variance(input[i]);
    }
  }

  void reflection_list_set_centroid_sq_width(af::ref<Reflection> r, 
      const af::const_ref< vec3<double> > &input) {
    DIALS_ASSERT(input.size() == r.size());        
    for (std::size_t i = 0; i < input.size(); ++i) {
      r[i].set_centroid_sq_width(input[i]);
    }
  }
  
  void reflection_list_set_intensity(af::ref<Reflection> r, 
      const af::const_ref<double> &input) {
    DIALS_ASSERT(input.size() == r.size());        
    for (std::size_t i = 0; i < input.size(); ++i) {
      r[i].set_intensity(input[i]);
    }
  }

  void reflection_list_set_intensity_variance(af::ref<Reflection> r, 
      const af::const_ref<double> &input) {
    DIALS_ASSERT(input.size() == r.size());      
    for (std::size_t i = 0; i < input.size(); ++i) {
      r[i].set_intensity_variance(input[i]);
    }
  }
  
  void reflection_list_set_corrected_intensity(af::ref<Reflection> r, 
      const af::const_ref<double> &input) {
    DIALS_ASSERT(input.size() == r.size());        
    for (std::size_t i = 0; i < input.size(); ++i) {
      r[i].set_corrected_intensity(input[i]);
    }
  }  
  
  void reflection_list_set_corrected_intensity_variance(af::ref<Reflection> r, 
      const af::const_ref<double> &input) {
    DIALS_ASSERT(input.size() == r.size());        
    for (std::size_t i = 0; i < input.size(); ++i) {
      r[i].set_corrected_intensity_variance(input[i]);
    }
  }  
    
  void reflection_list_set_shoebox(af::ref<Reflection> r, 
      const af::const_ref<double> &input) {
    std::size_t total_size = 0;
    for (std::size_t i = 0; i < r.size(); ++i) {
      int6 bbox = r[i].get_bounding_box();
      vec3<int> shape(bbox[5] - bbox[4], bbox[3] - bbox[2], bbox[1] - bbox[0]);
      DIALS_ASSERT(shape.const_ref().all_ge(0));
      total_size += shape[0] * shape[1] * shape[2];
    }
    DIALS_ASSERT(total_size == input.size());
    for (std::size_t i = 0, k = 0; i < r.size(); ++i) {
      int6 bbox = r[i].get_bounding_box();
      vec3<int> shape(bbox[5] - bbox[4], bbox[3] - bbox[2], bbox[1] - bbox[0]);    
      af::versa< double, af::c_grid<3> > s(af::c_grid<3>(shape[0], shape[1], shape[2]));
      for (std::size_t j = 0; j < s.size(); ++j) {
        s[j] = input[k++];
      }
      r[i].set_shoebox(s);
    }
  }  

  void reflection_list_set_shoebox_mask(af::ref<Reflection> r, 
      const af::const_ref<int> &input) {
    std::size_t total_size = 0;
    for (std::size_t i = 0; i < r.size(); ++i) {
      int6 bbox = r[i].get_bounding_box();
      vec3<int> shape(bbox[5] - bbox[4], bbox[3] - bbox[2], bbox[1] - bbox[0]);
      DIALS_ASSERT(shape.const_ref().all_ge(0));
      total_size += shape[0] * shape[1] * shape[2];
    }
    DIALS_ASSERT(total_size == input.size());
    for (std::size_t i = 0, k = 0; i < r.size(); ++i) {
      int6 bbox = r[i].get_bounding_box();
      vec3<int> shape(bbox[5] - bbox[4], bbox[3] - bbox[2], bbox[1] - bbox[0]);    
      af::versa< int, af::c_grid<3> > s(af::c_grid<3>(shape[0], shape[1], shape[2]));
      for (std::size_t j = 0; j < s.size(); ++j) {
        s[j] = input[k++];
      }
      r[i].set_shoebox_mask(s);
    }  
  }  
  
  void reflection_list_set_shoebox_background(af::ref<Reflection> r, 
      const af::const_ref<double> &input) {
    std::size_t total_size = 0;
    for (std::size_t i = 0; i < r.size(); ++i) {
      int6 bbox = r[i].get_bounding_box();
      vec3<int> shape(bbox[5] - bbox[4], bbox[3] - bbox[2], bbox[1] - bbox[0]);
      DIALS_ASSERT(shape.const_ref().all_ge(0));
      total_size += shape[0] * shape[1] * shape[2];
    }
    DIALS_ASSERT(total_size == input.size());
    for (std::size_t i = 0, k = 0; i < r.size(); ++i) {
      int6 bbox = r[i].get_bounding_box();
      vec3<int> shape(bbox[5] - bbox[4], bbox[3] - bbox[2], bbox[1] - bbox[0]);    
      af::versa< double, af::c_grid<3> > s(af::c_grid<3>(shape[0], shape[1], shape[2]));
      for (std::size_t j = 0; j < s.size(); ++j) {
        s[j] = input[k++];
      }
      r[i].set_shoebox_background(s);
    }   
  } 
   
  void reflection_list_set_transformed_shoebox(af::ref<Reflection> r, 
      const af::const_ref<double> &input) {
    if (input.size() == 0) {
      return;
    }
    DIALS_ASSERT(input.size() % r.size() == 0);
    std::size_t sz = (std::size_t)(std::exp(
      std::log(input.size() / r.size()) / 3.0) + 0.5);
    DIALS_ASSERT(r.size() * sz * sz * sz == input.size());    
    for (std::size_t i = 0, k = 0; i < r.size(); ++i) {
      af::versa< double, af::c_grid<3> > s(af::c_grid<3>(sz, sz, sz));
      for (std::size_t j = 0; j < s.size(); ++j) {
        s[j] = input[k++];
      }
      r[i].set_transformed_shoebox(s);
    }    
  }  
  
  static
  void set_shoebox(Reflection &obj, flex_double &data) {
    DIALS_ASSERT(data.accessor().all().size() == 3);
    obj.set_shoebox(af::versa<double, af::c_grid<3> >(
      data.handle(), af::c_grid<3>(data.accessor())));
  }

//  static
//  flex_double get_shoebox(Reflection &obj) {
//    return flex_double(obj.get_shoebox().handle(), 
//      obj.get_shoebox().accessor().as_flex_grid());
//  }
  
  static
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
  
  static
  void set_shoebox_background(Reflection &obj, flex_double &data) {
    DIALS_ASSERT(data.accessor().all().size() == 3);
    obj.set_shoebox_background(af::versa<double, af::c_grid<3> >(
      data.handle(), af::c_grid<3>(data.accessor())));
  }

//  static
//  flex_double get_shoebox_background(Reflection &obj) {
//    return flex_double(obj.get_shoebox_background().handle(), 
//      obj.get_shoebox_background().accessor().as_flex_grid());
//  }
  
  static
  void set_transformed_shoebox(Reflection &obj, flex_double &data) {
    DIALS_ASSERT(data.accessor().all().size() == 3);
    obj.set_transformed_shoebox(af::versa<double, af::c_grid<3> >(
      data.handle(), af::c_grid<3>(data.accessor())));
  }

//  static
//  flex_double get_transformed_shoebox(Reflection &obj) {
//    return flex_double(obj.get_transformed_shoebox().handle(), 
//      obj.get_transformed_shoebox().accessor().as_flex_grid());
//  }
  
  static
  void set_transformed_shoebox_background(Reflection &obj, flex_double &data) {
    DIALS_ASSERT(data.accessor().all().size() == 3);
    obj.set_transformed_shoebox_background(af::versa<double, af::c_grid<3> >(
      data.handle(), af::c_grid<3>(data.accessor())));
  }

//  static
//  flex_double get_transformed_shoebox_background(Reflection &obj) {
//    return flex_double(obj.get_transformed_shoebox_background().handle(), 
//      obj.get_transformed_shoebox_background().accessor().as_flex_grid());
//  }
  
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
      .def("__str__", &reflection_to_string)
      .def_pickle(reflection::ReflectionPickleSuite());          

    scitbx::af::boost_python::flex_wrapper 
        <int6>::plain("flex_int6");

    scitbx::af::boost_python::flex_wrapper 
      <Reflection, return_internal_reference<> >::plain("ReflectionList")
        .def("__init__", make_constructor(&init_from_observation_and_shoebox))
        .def_pickle(scitbx::af::boost_python::flex_pickle_double_buffered<
          Reflection, reflection::to_string, reflection::from_string>())      
        .def("miller_index", &reflection_list_get_miller_index)
        .def("miller_index", &reflection_list_set_miller_index)
        .def("status", &reflection_list_get_status)
        .def("status", &reflection_list_set_status)
        .def("entering", &reflection_list_get_entering)
        .def("entering", &reflection_list_set_entering)
        .def("rotation_angle", &reflection_list_get_rotation_angle)
        .def("rotation_angle", &reflection_list_set_rotation_angle)
        .def("beam_vector", &reflection_list_get_beam_vector)
        .def("beam_vector", &reflection_list_set_beam_vector)
        .def("image_coord_mm", &reflection_list_get_image_coord_mm)
        .def("image_coord_mm", &reflection_list_set_image_coord_mm)
        .def("image_coord_px", &reflection_list_get_image_coord_px)
        .def("image_coord_px", &reflection_list_set_image_coord_px)
        .def("frame_number", &reflection_list_get_frame_number)
        .def("frame_number", &reflection_list_set_frame_number)
        .def("panel_number", &reflection_list_get_panel_number)
        .def("panel_number", &reflection_list_set_panel_number)
        .def("bounding_box", &reflection_list_get_bounding_box)
        .def("bounding_box", &reflection_list_set_bounding_box)
        .def("centroid_position", &reflection_list_get_centroid_position)
        .def("centroid_position", &reflection_list_set_centroid_position)
        .def("centroid_variance", &reflection_list_get_centroid_variance)
        .def("centroid_variance", &reflection_list_set_centroid_variance)
        .def("centroid_sq_width", &reflection_list_get_centroid_sq_width)
        .def("centroid_sq_width", &reflection_list_set_centroid_sq_width)
        .def("intensity", &reflection_list_get_intensity)
        .def("intensity", &reflection_list_set_intensity)
        .def("intensity_variance", &reflection_list_get_intensity_variance)
        .def("intensity_variance", &reflection_list_set_intensity_variance)
        .def("corrected_intensity", &reflection_list_get_corrected_intensity)
        .def("corrected_intensity", &reflection_list_set_corrected_intensity)        
        .def("corrected_intensity_variance", 
          &reflection_list_get_corrected_intensity_variance)
        .def("corrected_intensity_variance", 
          &reflection_list_set_corrected_intensity_variance)     
        .def("shoebox", &reflection_list_get_shoebox)
        .def("shoebox", &reflection_list_set_shoebox)
        .def("shoebox_mask", &reflection_list_get_shoebox_mask)
        .def("shoebox_mask", &reflection_list_set_shoebox_mask)
        .def("shoebox_background", &reflection_list_get_shoebox_background)
        .def("shoebox_background", &reflection_list_set_shoebox_background)
        .def("transformed_shoebox", &reflection_list_get_transformed_shoebox)
        .def("transformed_shoebox", &reflection_list_set_transformed_shoebox);
  }

}}} // namespace dials::model::boost_python

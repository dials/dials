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
//#include <scitbx/serialization/single_buffered.h>
#include <scitbx/array_family/ref_reductions.h>
#include <dials/model/data/reflection.h>
#include <dials/error.h>
//#include <scitbx/array_family/boost_python/ref_pickle_double_buffered.h>
#include <dials/model/data/boost_python/pickle.h>
#include <scitbx/array_family/boost_python/flex_pickle_double_buffered.h>




//#include <scitbx/array_family/boost_python/flex_pickle_single_buffered.h>

namespace dials { namespace model { namespace boost_python {

  //using namespace scitbx::af::boost_python;
 // using namespace scitbx::serialization::single_buffered;


  using namespace boost::python;
 // using scitbx::af::boost_python::flex_pickle_single_buffered;


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
      r.set_shoebox(extract<const flex_double&>(state[12]));
      r.set_shoebox_mask(extract<const flex_int&>(state[13]));
      r.set_transformed_shoebox(extract<const flex_double&>(state[14]));
    }
  };

  using scitbx::af::flex_grid;
  using scitbx::af::flex_bool;
//  using scitbx::af::boost_python::flex_pickle_single_buffered;
//  using scitbx::af::boost_python::pickle_size_per_element;
  typedef scitbx::af::flex<cctbx::miller::index<> >::type flex_miller_index;
  typedef scitbx::af::flex<vec2<double> >::type flex_vec2_double;
  typedef scitbx::af::flex<vec3<double> >::type flex_vec3_double;
  typedef scitbx::af::flex<vec3<int> >::type flex_vec3_int;
  typedef scitbx::af::flex<int6>::type flex_int6;

  flex_miller_index reflection_list_get_miller_index(const ReflectionList &r) {
    flex_miller_index result(r.size());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = r[i].get_miller_index();
    }
    return result;
  }
  
  flex_int reflection_list_get_status(const ReflectionList &r) {
    flex_int result(r.size());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = r[i].get_status();
    }
    return result;
  }
  
  flex_bool reflection_list_get_entering(const ReflectionList &r) {
    flex_bool result(r.size());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = r[i].get_entering();
    }
    return result;
  }

  flex_double reflection_list_get_rotation_angle(const ReflectionList &r) {
    flex_double result(r.size());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = r[i].get_rotation_angle();
    }
    return result;
  }

  flex_vec3_double reflection_list_get_beam_vector(const ReflectionList &r) {
    flex_vec3_double result(r.size());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = r[i].get_beam_vector();
    }
    return result;
  }
  
  flex_vec2_double reflection_list_get_image_coord_mm(const ReflectionList &r) {
    flex_vec2_double result(r.size());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = r[i].get_image_coord_mm();
    }
    return result;
  }
  
  flex_vec2_double reflection_list_get_image_coord_px(const ReflectionList &r) {
    flex_vec2_double result(r.size());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = r[i].get_image_coord_px();
    }
    return result;
  }
  
  flex_double reflection_list_get_frame_number(const ReflectionList &r) {
    flex_double result(r.size());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = r[i].get_frame_number();
    }
    return result;
  }
  
  flex_int reflection_list_get_panel_number(const ReflectionList &r) {
    flex_int result(r.size());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = r[i].get_panel_number();
    }
    return result;
  }
  
  flex_int6 reflection_list_get_bounding_box(const ReflectionList &r) {
    flex_int6 result(r.size());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = r[i].get_bounding_box();
    }
    return result;
  }
  
  flex_vec3_double reflection_list_get_centroid_position(const ReflectionList &r) {
    flex_vec3_double result(r.size());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = r[i].get_centroid_position();
    }
    return result;
  }
  
  flex_vec3_double reflection_list_get_centroid_variance(const ReflectionList &r) {
    flex_vec3_double result(r.size());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = r[i].get_centroid_variance();
    }
    return result;
  }
  
  flex_vec3_double reflection_list_get_centroid_sq_width(const ReflectionList &r) {
    flex_vec3_double result(r.size());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = r[i].get_centroid_sq_width();
    }
    return result;
  }
  
  flex_double reflection_list_get_intensity(const ReflectionList &r) {
    flex_double result(r.size());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = r[i].get_intensity();
    }
    return result;
  }
  
  flex_double reflection_list_get_intensity_variance(const ReflectionList &r) {
    flex_double result(r.size());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = r[i].get_intensity_variance();
    }
    return result;
  }  
  
 
  flex_double reflection_list_get_shoebox(const ReflectionList &r) {
    std::size_t result_size = 0;
    for (std::size_t i = 0; i < r.size(); ++i) {
      result_size += r[i].get_shoebox().size();
    }
    flex_double result(result_size);
    for (std::size_t i = 0, k = 0; i < r.size(); ++i) {
      flex_double s = r[i].get_shoebox();
      for (std::size_t j = 0; j < s.size(); ++j) {
        result[k++] = s[j];
      }
    }
    return result;
  }
  
  flex_int reflection_list_get_shoebox_mask(const ReflectionList &r) {
    std::size_t result_size = 0;
    for (std::size_t i = 0; i < r.size(); ++i) {
      result_size += r[i].get_shoebox_mask().size();
    }
    flex_int result(result_size);
    for (std::size_t i = 0, k = 0; i < r.size(); ++i) {
      flex_int s = r[i].get_shoebox_mask();
      for (std::size_t j = 0; j < s.size(); ++j) {
        result[k++] = s[j];
      }
    }
    return result;
  }
  
  flex_double reflection_list_get_shoebox_background(const ReflectionList &r) {
    std::size_t result_size = 0;
    for (std::size_t i = 0; i < r.size(); ++i) {
      result_size += r[i].get_shoebox_background().size();
    }
    flex_double result(result_size);
    for (std::size_t i = 0, k = 0; i < r.size(); ++i) {
      flex_double s = r[i].get_shoebox_background();
      for (std::size_t j = 0; j < s.size(); ++j) {
        result[k++] = s[j];
      }
    }
    return result;
  }
  
  flex_double reflection_list_get_transformed_shoebox(const ReflectionList &r) {
    std::size_t result_size = 0;
    for (std::size_t i = 0; i < r.size(); ++i) {
      result_size += r[i].get_transformed_shoebox().size();
    }
    flex_double result(result_size);
    for (std::size_t i = 0, k = 0; i < r.size(); ++i) {
      flex_double s = r[i].get_transformed_shoebox();
      for (std::size_t j = 0; j < s.size(); ++j) {
        result[k++] = s[j];
      }
    }
    return result;
  }
  
  
  void reflection_list_set_miller_index(ReflectionList &r, 
      const flex_miller_index &input) {
    DIALS_ASSERT(input.size() == r.size());        
    for (std::size_t i = 0; i < input.size(); ++i) {
      r[i].set_miller_index(input[i]);
    }
  }

  void reflection_list_set_status(ReflectionList &r, 
      const flex_int &input) {
    DIALS_ASSERT(input.size() == r.size());        
    for (std::size_t i = 0; i < input.size(); ++i) {
      r[i].set_status(input[i]);
    }
  }

  void reflection_list_set_entering(ReflectionList &r, 
      const flex_bool &input) {
    DIALS_ASSERT(input.size() == r.size());        
    for (std::size_t i = 0; i < input.size(); ++i) {
      r[i].set_entering(input[i]);
    }
  }


  void reflection_list_set_rotation_angle(ReflectionList &r, 
      const flex_double &input) {
    DIALS_ASSERT(input.size() == r.size());        
    for (std::size_t i = 0; i < input.size(); ++i) {
      r[i].set_rotation_angle(input[i]);
    }
  }

  void reflection_list_set_beam_vector(ReflectionList &r, 
      const flex_vec3_double &input) {
    DIALS_ASSERT(input.size() == r.size());        
    for (std::size_t i = 0; i < input.size(); ++i) {
      r[i].set_beam_vector(input[i]);
    }
  }

  void reflection_list_set_image_coord_mm(ReflectionList &r, 
      const flex_vec2_double &input) {
    DIALS_ASSERT(input.size() == r.size());        
    for (std::size_t i = 0; i < input.size(); ++i) {
      r[i].set_image_coord_mm(input[i]);
    }
  }

  void reflection_list_set_image_coord_px(ReflectionList &r, 
      const flex_vec2_double &input) {
    DIALS_ASSERT(input.size() == r.size());        
    for (std::size_t i = 0; i < input.size(); ++i) {
      r[i].set_image_coord_px(input[i]);
    }
  }
  
  void reflection_list_set_frame_number(ReflectionList &r, 
      const flex_double &input) {
    DIALS_ASSERT(input.size() == r.size());        
    for (std::size_t i = 0; i < input.size(); ++i) {
      r[i].set_frame_number(input[i]);
    }
  }

  void reflection_list_set_panel_number(ReflectionList &r, 
      const flex_int &input) {
    DIALS_ASSERT(input.size() == r.size());        
    for (std::size_t i = 0; i < input.size(); ++i) {
      r[i].set_panel_number(input[i]);
    }
  }

  void reflection_list_set_bounding_box(ReflectionList &r, 
      const flex_int6 &input) {
    DIALS_ASSERT(input.size() == r.size());        
    for (std::size_t i = 0; i < input.size(); ++i) {
      r[i].set_bounding_box(input[i]);
    }
  }
  
  void reflection_list_set_centroid_position(ReflectionList &r, 
      const flex_vec3_double &input) {
    DIALS_ASSERT(input.size() == r.size());        
    for (std::size_t i = 0; i < input.size(); ++i) {
      r[i].set_centroid_position(input[i]);
    }
  }
  
  void reflection_list_set_centroid_variance(ReflectionList &r, 
      const flex_vec3_double &input) {
    DIALS_ASSERT(input.size() == r.size());        
    for (std::size_t i = 0; i < input.size(); ++i) {
      r[i].set_centroid_variance(input[i]);
    }
  }

  void reflection_list_set_centroid_sq_width(ReflectionList &r, 
      const flex_vec3_double &input) {
    DIALS_ASSERT(input.size() == r.size());        
    for (std::size_t i = 0; i < input.size(); ++i) {
      r[i].set_centroid_sq_width(input[i]);
    }
  }
  
  void reflection_list_set_intensity(ReflectionList &r, 
      const flex_double &input) {
    DIALS_ASSERT(input.size() == r.size());        
    for (std::size_t i = 0; i < input.size(); ++i) {
      r[i].set_intensity(input[i]);
    }
  }

  void reflection_list_set_intensity_variance(ReflectionList &r, 
      const flex_double &input) {
    DIALS_ASSERT(input.size() == r.size());      
    for (std::size_t i = 0; i < input.size(); ++i) {
      r[i].set_intensity_variance(input[i]);
    }
  }
  
  void reflection_list_set_shoebox(ReflectionList &r, 
      const flex_double &input) {
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
      flex_double s(flex_grid<>(shape[0], shape[1], shape[2]));
      for (std::size_t j = 0; j < s.size(); ++j) {
        s[j] = input[k++];
      }
      r[i].set_shoebox(s);
    }
  }  

  void reflection_list_set_shoebox_mask(ReflectionList &r, 
      const flex_int &input) {
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
      flex_int s(flex_grid<>(shape[0], shape[1], shape[2]));
      for (std::size_t j = 0; j < s.size(); ++j) {
        s[j] = input[k++];
      }
      r[i].set_shoebox_mask(s);
    }  
  }  
  
  void reflection_list_set_shoebox_background(ReflectionList &r, 
      const flex_double &input) {
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
      flex_double s(flex_grid<>(shape[0], shape[1], shape[2]));
      for (std::size_t j = 0; j < s.size(); ++j) {
        s[j] = input[k++];
      }
      r[i].set_shoebox_background(s);
    }   
  } 
   
  void reflection_list_set_transformed_shoebox(ReflectionList &r, 
      const flex_double &input) {
    if (input.size() == 0) {
      return;
    }
    DIALS_ASSERT(input.size() % r.size() == 0);
    std::size_t sz = (std::size_t)(std::exp(
      std::log(input.size() / r.size()) / 3.0) + 0.5);
    DIALS_ASSERT(r.size() * sz * sz * sz == input.size());    
    for (std::size_t i = 0, k = 0; i < r.size(); ++i) {
      flex_double s(flex_grid<>(sz, sz, sz));
      for (std::size_t j = 0; j < s.size(); ++j) {
        s[j] = input[k++];
      }
      r[i].set_transformed_shoebox(s);
    }    
  }  
  
  
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
      .add_property("shoebox",
        &Reflection::get_shoebox,
        &Reflection::set_shoebox)
      .add_property("shoebox_mask",
        &Reflection::get_shoebox_mask,
        &Reflection::set_shoebox_mask)
      .add_property("shoebox_background",
        &Reflection::get_shoebox_background,
        &Reflection::set_shoebox_background)        
      .add_property("transformed_shoebox",
        &Reflection::get_transformed_shoebox,
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

//    scitbx::af::boost_python::flex_wrapper 
//        <int6>::plain("flex_int6")
//      .def_pickle(flex_pickle_single_buffered<int6,
//        6*pickle_size_per_element<int>::value>());

    scitbx::af::boost_python::flex_wrapper 
      <Reflection, return_internal_reference<> >::plain("ReflectionList")
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

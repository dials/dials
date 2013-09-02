/*
 * prediction.cc
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
#include <dials/model/data/prediction.h>

namespace dials { namespace model { namespace boost_python {

  using namespace boost::python;
 
  static
  vec3<double> position_data_get_px(const Prediction::PositionData& obj) {
    return obj.px;
  }

  static
  void position_data_set_px(Prediction::PositionData& obj, vec3<double> v) {
    obj.px = v;
  }

  static
  vec3<double> position_data_get_mm(const Prediction::PositionData& obj) {
    return obj.mm;
  }
  
  static
  void position_data_set_mm(Prediction::PositionData& obj, vec3<double> v) {
    obj.mm = v;
  }
 
  static
  MillerIndex prediction_get_miller_index(const Prediction &obj) {
    return obj.miller_index;
  }
 
  static
  void prediction_set_miller_index(Prediction &obj, MillerIndex v) {
    obj.miller_index = v;
  }
  
  static
  vec3<double> prediction_get_beam_vector(const Prediction& obj) {
    return obj.beam_vector;
  }
  
  static
  void prediction_set_beam_vector(Prediction& obj, vec3<double> v) {
    obj.beam_vector = v;
  }

  void export_prediction()
  {
    class_<Prediction::PositionData>("PositionData")
      .add_property("px",
        &position_data_get_px,
        &position_data_set_px) 
      .add_property("mm",
        &position_data_get_mm,
        &position_data_set_mm);     
    
    class_<Prediction>("Prediction")
      .add_property("miller_index",
        &prediction_get_miller_index,
        &prediction_set_miller_index)
      .add_property("beam_vector",
        &prediction_get_beam_vector,
        &prediction_set_beam_vector)        
      .def_readwrite("position", &Prediction::position)
      .def_readwrite("panel", &Prediction::panel)
      .def_readwrite("entering", &Prediction::entering);
  }

}}} // namespace dials::model::boost_python

/*
 * observation.cc
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
#include <dials/model/data/observation.h>

namespace dials { namespace model { namespace boost_python {

  using namespace boost::python;

  static
  vec3<double> position_data_get_position(const Position::PositionData &obj) {
    return obj.position;
  }

  static
  vec3<double> position_data_get_variance(const Position::PositionData &obj) {
    return obj.variance;
  }
  
  static
  vec3<double> position_data_get_std_err_sq(const Position::PositionData &obj) {
    return obj.std_err_sq;
  }

  static
  void position_data_set_position(Position::PositionData &obj, vec3<double> v) {
    obj.position = v;
  }

  static
  void position_data_set_variance(Position::PositionData &obj, vec3<double> v) {
    obj.variance = v;
  }
  
  static
  void position_data_set_std_err_sq(Position::PositionData &obj, vec3<double> v) {
    obj.std_err_sq = v;
  }
  
  void export_observation()
  {
    class_<Intensity::IntensityData>("IntensityData")
      .def(init<double, double>((
        arg("value"),
        arg("variance"))))
      .def_readwrite("value", &Intensity::IntensityData::value)
      .def_readwrite("variance", &Intensity::IntensityData::variance);
    
    class_<Intensity>("Intensity")
      .def(init<double, double>((
        arg("observed_value"),
        arg("observed_variance"))))
      .def(init<double, double, double, double>((
        arg("observed_value"),
        arg("observed_variance"),
        arg("corrected_value"),
        arg("corrected_variance"))))
      .def(init<const Intensity::IntensityData&>((
        arg("observed"))))
      .def(init<
          const Intensity::IntensityData&, 
          const Intensity::IntensityData&>((
        arg("observed"),
        arg("corrected"))))
      .def_readwrite("observed", &Intensity::observed)
      .def_readwrite("corrected", &Intensity::corrected);

    class_<Position::PositionData>("PositionData")
      .def(init<vec3<double>, vec3<double>, vec3<double> >((
        arg("position"),
        arg("variance"),
        arg("std_err_sq"))))
      .add_property("position", 
        &position_data_get_position,
        &position_data_set_position)
      .add_property("variance", 
        &position_data_get_variance,
        &position_data_set_variance)
      .add_property("std_err_sq",
        &position_data_get_std_err_sq,
        &position_data_set_std_err_sq);
      
    class_<Position>("Position")
      .def(init<vec3<double>, vec3<double>, vec3<double> >((
        arg("px_position"),
        arg("px_variance"),
        arg("px_std_err_sq"))))
      .def(init<vec3<double>, vec3<double>, vec3<double>,
                vec3<double>, vec3<double>, vec3<double> >((
        arg("px_position"),
        arg("px_variance"),
        arg("px_std_err_sq"),
        arg("mm_position"),
        arg("mm_variance"),
        arg("mm_std_err_sq"))))   
      .def(init<const Position::PositionData&>((
        arg("px"))))
      .def(init<const Position::PositionData&, const Position::PositionData&>((
        arg("px"),
        arg("mm"))))    
      .def_readwrite("px", &Position::px)
      .def_readwrite("mm", &Position::mm);
  
    class_<Observation>("Observation")
      .def(init<const Position&>((
        arg("centroid"))))
      .def(init<const Intensity&>((
        arg("intensity"))))
      .def(init<const Position&, const Intensity&>((
        arg("centroid"),
        arg("intensity"))))
      .def_readwrite("centroid", &Observation::centroid)
      .def_readwrite("intensity", &Observation::intensity);
  }

}}} // namespace dials::model::boost_python

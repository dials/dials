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

  /** Wrapper function to get the millimeter x/y coordinate */
  static vec2<double> position_data_get_mm_xy(const Prediction::PositionData& obj) {
    return vec2<double>(obj.mm[0], obj.mm[1]);
  }

  /** Wrapper function to set the millimeter x/y coordinate */
  static void position_data_set_mm_xy(Prediction::PositionData& obj, vec2<double> v) {
    obj.mm[0] = v[0];
    obj.mm[1] = v[1];
  }

  /** Wrapper function to get the angle */
  static double position_data_get_angle(const Prediction::PositionData& obj) {
    return obj.mm[2];
  }

  /** Wrapper function to set the angle */
  static void position_data_set_angle(Prediction::PositionData& obj, double v) {
    obj.mm[2] = v;
  }

  /** Wrapper function to get the pixel x/y coordinate */
  static vec2<double> position_data_get_px_xy(const Prediction::PositionData& obj) {
    return vec2<double>(obj.px[0], obj.px[1]);
  }

  /** Wrapper function to set the x/y coordinate */
  static void position_data_set_px_xy(Prediction::PositionData& obj, vec2<double> v) {
    obj.px[0] = v[0];
    obj.px[1] = v[1];
  }

  /** Wrapper function to get the frame */
  static double position_data_get_frame(const Prediction::PositionData& obj) {
    return obj.px[2];
  }

  /** Wrapper function to set the frame */
  static void position_data_set_frame(Prediction::PositionData& obj, double v) {
    obj.px[2] = v;
  }

  void export_prediction() {
    class_<Prediction::PositionData>("PositionData")
      .add_property("px",
                    make_getter(&Prediction::PositionData::px,
                                return_value_policy<return_by_value>()),
                    make_setter(&Prediction::PositionData::px,
                                return_value_policy<return_by_value>()))
      .add_property("mm",
                    make_getter(&Prediction::PositionData::mm,
                                return_value_policy<return_by_value>()),
                    make_setter(&Prediction::PositionData::mm,
                                return_value_policy<return_by_value>()))
      .add_property("mm_xy", &position_data_get_mm_xy, &position_data_set_mm_xy)
      .add_property("angle", &position_data_get_angle, &position_data_set_angle)
      .add_property("px_xy", &position_data_get_px_xy, &position_data_set_px_xy)
      .add_property("frame", &position_data_get_frame, &position_data_set_frame)
      .def("__eq__", &Prediction::PositionData::operator==)
      .def("__ne__", &Prediction::PositionData::operator!=);

    class_<Prediction>("Prediction")
      .def(init<const Prediction&>())
      .add_property(
        "miller_index",
        make_getter(&Prediction::miller_index, return_value_policy<return_by_value>()),
        make_setter(&Prediction::miller_index, return_value_policy<return_by_value>()))
      .add_property(
        "beam_vector",
        make_getter(&Prediction::beam_vector, return_value_policy<return_by_value>()),
        make_setter(&Prediction::beam_vector, return_value_policy<return_by_value>()))
      .def_readwrite("position", &Prediction::position)
      .def_readwrite("panel", &Prediction::panel)
      .def_readwrite("entering", &Prediction::entering)
      .def_readwrite("crystal", &Prediction::crystal)
      .def("__eq__", &Prediction::operator==)
      .def("__ne__", &Prediction::operator!=);
  }

}}}  // namespace dials::model::boost_python

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
#include <scitbx/vec2.h>
#include <scitbx/vec3.h>
#include <dials/model/data/observation.h>

namespace dials { namespace model { namespace boost_python {

  using namespace boost::python;
  using scitbx::vec2;
  using scitbx::vec3;

  /** Wrapper function to get the xy pixel coordinate */
  static vec2<double> centroid_get_px_xy(const Centroid &obj) {
    return vec2<double>(obj.px.position[0], obj.px.position[1]);
  }

  /** Wrapper function to set the xy pixel coordinate */
  static void centroid_set_px_xy(Centroid &obj, vec2<double> v) {
    obj.px.position[0] = v[0];
    obj.px.position[1] = v[1];
  }

  /** Wrapper function to get the frame */
  static double centroid_get_frame(const Centroid &obj) {
    return obj.px.position[2];
  }

  /** Wrapper function to set the frame */
  static void centroid_set_frame(Centroid &obj, double v) {
    obj.px.position[2] = v;
  }

  /** Wrapper function to get the millimeter xy coordinate */
  static vec2<double> centroid_get_mm_xy(const Centroid &obj) {
    return vec2<double>(obj.mm.position[0], obj.mm.position[1]);
  }

  /** Wrapper function to set the millimeter xy coordinate */
  static void centroid_set_mm_xy(Centroid &obj, vec2<double> v) {
    obj.mm.position[0] = v[0];
    obj.mm.position[1] = v[1];
  }

  /** Wrapper function to get the angle */
  static double centroid_get_angle(const Centroid &obj) {
    return obj.mm.position[2];
  }

  /** Wrapper function to set the angle */
  static void centroid_set_angle(Centroid &obj, double v) {
    obj.mm.position[2] = v;
  }

  void export_observation() {
    class_<Intensity::IntensityData>("IntensityData")
      .def(init<double, double, bool>((arg("value"), arg("variance"), arg("success"))))
      .def_readwrite("value", &Intensity::IntensityData::value)
      .def_readwrite("variance", &Intensity::IntensityData::variance)
      .def_readwrite("success", &Intensity::IntensityData::success)
      .def("__eq__", &Intensity::IntensityData::operator==)
      .def("__ne__", &Intensity::IntensityData::operator!=);

    class_<Intensity>("Intensity")
      .def(init<double, double, bool>(
        (arg("observed_value"), arg("observed_variance"), arg("observed_success"))))
      .def(init<double, double, bool, double, double, bool>((arg("observed_value"),
                                                             arg("observed_variance"),
                                                             arg("observed_success"),
                                                             arg("corrected_value"),
                                                             arg("corrected_variance"),
                                                             arg("corrected_success"))))
      .def(init<const Intensity::IntensityData &>((arg("observed"))))
      .def(init<const Intensity::IntensityData &, const Intensity::IntensityData &>(
        (arg("observed"), arg("corrected"))))
      .def_readwrite("observed", &Intensity::observed)
      .def_readwrite("corrected", &Intensity::corrected)
      .def_readwrite("background", &Intensity::background)
      .def("__eq__", &Intensity::operator==)
      .def("__ne__", &Intensity::operator!=);

    class_<Centroid::CentroidData>("CentroidData")
      .def(init<vec3<double>, vec3<double>, vec3<double> >(
        (arg("position"), arg("variance"), arg("std_err_sq"))))
      .add_property("position",
                    make_getter(&Centroid::CentroidData::position,
                                return_value_policy<return_by_value>()),
                    make_setter(&Centroid::CentroidData::position,
                                return_value_policy<return_by_value>()))
      .add_property("variance",
                    make_getter(&Centroid::CentroidData::variance,
                                return_value_policy<return_by_value>()),
                    make_setter(&Centroid::CentroidData::variance,
                                return_value_policy<return_by_value>()))
      .add_property("std_err_sq",
                    make_getter(&Centroid::CentroidData::std_err_sq,
                                return_value_policy<return_by_value>()),
                    make_setter(&Centroid::CentroidData::std_err_sq,
                                return_value_policy<return_by_value>()))
      .def("__eq__", &Centroid::CentroidData::operator==)
      .def("__ne__", &Centroid::CentroidData::operator!=);

    class_<Centroid>("Centroid")
      .def(init<vec3<double>, vec3<double>, vec3<double> >(
        (arg("px_position"), arg("px_variance"), arg("px_std_err_sq"))))
      .def(init<vec3<double>,
                vec3<double>,
                vec3<double>,
                vec3<double>,
                vec3<double>,
                vec3<double> >((arg("px_position"),
                                arg("px_variance"),
                                arg("px_std_err_sq"),
                                arg("mm_position"),
                                arg("mm_variance"),
                                arg("mm_std_err_sq"))))
      .def(init<const Centroid::CentroidData &>((arg("px"))))
      .def(init<const Centroid::CentroidData &, const Centroid::CentroidData &>(
        (arg("px"), arg("mm"))))
      .def_readwrite("px", &Centroid::px)
      .def_readwrite("mm", &Centroid::mm)
      .add_property("px_xy", &centroid_get_px_xy, &centroid_set_px_xy)
      .add_property("frame", &centroid_get_frame, &centroid_set_frame)
      .add_property("mm_xy", &centroid_get_mm_xy, &centroid_set_mm_xy)
      .add_property("angle", &centroid_get_angle, &centroid_set_angle)
      .def("update_mm",
           (void (Centroid::*)(const Detector &, const Scan &)) & Centroid::update_mm,
           (arg("detector"), arg("scan")))
      .def("update_mm",
           (void (Centroid::*)(std::size_t, const Detector &, const Scan &))
             & Centroid::update_mm,
           (arg("panel"), arg("detector"), arg("scan")))
      .def("resolution",
           (double (Centroid::*)(const BeamBase &, const Detector &) const)
             & Centroid::resolution,
           (arg("beam"), arg("detector")))
      .def("resolution",
           (double (Centroid::*)(std::size_t, const BeamBase &, const Detector &) const)
             & Centroid::resolution,
           (arg("panel"), arg("beam"), arg("detector")))
      .def("__eq__", &Centroid::operator==)
      .def("__ne__", &Centroid::operator!=);

    class_<Observation>("Observation")
      .def(init<const Observation &>())
      .def(init<const Centroid &>((arg("centroid"))))
      .def(init<const Intensity &>((arg("intensity"))))
      .def(
        init<const Centroid &, const Intensity &>((arg("centroid"), arg("intensity"))))
      .def(init<std::size_t>((arg("panel"))))
      .def(init<std::size_t, const Centroid &>((arg("panel"), arg("centroid"))))
      .def(init<std::size_t, const Intensity &>((arg("panel"), arg("intensity"))))
      .def(init<std::size_t, const Centroid &, const Intensity &>(
        (arg("panel"), arg("centroid"), arg("intensity"))))
      .def_readwrite("panel", &Observation::panel)
      .def_readwrite("centroid", &Observation::centroid)
      .def_readwrite("intensity", &Observation::intensity)
      .def("update_centroid_mm",
           &Observation::update_centroid_mm,
           (arg("detector"), arg("scan")))
      .def("resolution", &Observation::resolution, (arg("beam"), arg("detector")))
      .def("__eq__", &Observation::operator==)
      .def("__ne__", &Observation::operator!=);
  }

}}}  // namespace dials::model::boost_python

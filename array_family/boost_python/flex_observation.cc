/*
 * flex_observation.cc
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
#include <scitbx/array_family/boost_python/flex_wrapper.h>
#include <scitbx/array_family/ref_reductions.h>
#include <dials/model/data/observation.h>

namespace dials { namespace af { namespace boost_python {

  using namespace boost::python;
  
  using scitbx::af::ref;
  using scitbx::af::const_ref;
  using scitbx::af::shared;
  using scitbx::vec2;
  using scitbx::vec3;
  using dxtbx::model::Detector;
  using dxtbx::model::Scan;
  using dials::model::Observation;

  static
  shared<std::size_t> observation_panel(
      const const_ref<Observation> &obj) {
     shared<std::size_t> result(obj.size());
     for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = obj[i].panel;
     }
     return result;
  }

  static
  shared< vec3<double> > observation_centroid_px_position(
      const const_ref<Observation> &obj) {
     shared< vec3<double> > result(obj.size());
     for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = obj[i].centroid.px.position;
     }
     return result;
  }

  static
  shared< vec3<double> > observation_centroid_px_variance(
      const const_ref<Observation> &obj) {
     shared< vec3<double> > result(obj.size());
     for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = obj[i].centroid.px.variance;
     }
     return result;
  }

  static
  shared< vec3<double> > observation_centroid_px_std_err_sq(
      const const_ref<Observation> &obj) {
     shared< vec3<double> > result(obj.size());
     for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = obj[i].centroid.px.std_err_sq;
     }
     return result;
  }

  static
  shared< vec3<double> > observation_centroid_mm_position(
      const const_ref<Observation> &obj) {
     shared< vec3<double> > result(obj.size());
     for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = obj[i].centroid.mm.position;
     }
     return result;
  }

  static
  shared< vec3<double> > observation_centroid_mm_variance(
      const const_ref<Observation> &obj) {
     shared< vec3<double> > result(obj.size());
     for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = obj[i].centroid.mm.variance;
     }
     return result;
  }

  static
  shared< vec3<double> > observation_centroid_mm_std_err_sq(
      const const_ref<Observation> &obj) {
     shared< vec3<double> > result(obj.size());
     for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = obj[i].centroid.mm.std_err_sq;
     }
     return result;
  }

  static
  shared< vec2<double> > observation_centroid_px_position_xy(
      const const_ref<Observation> &obj) {
     shared<vec2<double> > result(obj.size());
     for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = vec2<double>(obj[i].centroid.px.position[0], 
                               obj[i].centroid.px.position[1]);
     }
     return result;
  }

  static
  shared<double> observation_centroid_position_frame(
      const const_ref<Observation> &obj) {
     shared<double> result(obj.size());
     for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = obj[i].centroid.px.position[2];
     }
     return result;
  }

  static
  shared< vec2<double> > observation_centroid_mm_position_xy(
      const const_ref<Observation> &obj) {
     shared<vec2<double> > result(obj.size());
     for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = vec2<double>(obj[i].centroid.mm.position[0], 
                               obj[i].centroid.mm.position[1]);
     }
     return result;
  }

  static
  shared<double> observation_centroid_position_angle(
      const const_ref<Observation> &obj) {
     shared<double> result(obj.size());
     for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = obj[i].centroid.mm.position[2];
     }
     return result;
  }

  static
  shared<double> observation_intensity_observed_value(
      const const_ref<Observation> &obj) {
     shared<double> result(obj.size());
     for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = obj[i].intensity.observed.value;
     }
     return result;
  }

  static
  shared<double> observation_intensity_observed_variance(
      const const_ref<Observation> &obj) {
     shared<double> result(obj.size());
     for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = obj[i].intensity.observed.variance;
     }
     return result;
  }

  static
  shared<double> observation_intensity_corrected_value(
      const const_ref<Observation> &obj) {
     shared<double> result(obj.size());
     for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = obj[i].intensity.corrected.value;
     }
     return result;
  }

  static
  shared<double> observation_intensity_corrected_variance(
      const const_ref<Observation> &obj) {
     shared<double> result(obj.size());
     for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = obj[i].intensity.corrected.variance;
     }
     return result;
  }
  
  void observation_update_centroid_mm(ref<Observation> &obj, 
      const Detector &d, const Scan &s) {
    for (std::size_t i = 0; i < obj.size(); ++i) {
      obj[i].update_centroid_mm(d, s);
    }
  }
  
  void export_flex_observation()
  {
    scitbx::af::boost_python::flex_wrapper <
        Observation, return_internal_reference<> >::plain("observation")
      .def("panel",
        &observation_panel)
      .def("centroid_px_position",
        &observation_centroid_px_position)
      .def("centroid_px_variance",
        &observation_centroid_px_variance)
      .def("centroid_px_std_err_eq",
        &observation_centroid_px_std_err_sq)
      .def("centroid_mm_position",
        &observation_centroid_mm_position)
      .def("centroid_mm_variance",
        &observation_centroid_mm_variance)
      .def("centroid_mm_std_err_eq",
        &observation_centroid_mm_std_err_sq)
      .def("centroid_px_position_xy",
        &observation_centroid_px_position_xy)
      .def("centroid_position_frame",
        &observation_centroid_position_frame)
      .def("centroid_mm_position_xy",
        &observation_centroid_mm_position_xy)
      .def("centroid_position_angle",
        &observation_centroid_position_angle)
      .def("intensity_observed_value",
        &observation_intensity_observed_value)
      .def("intensity_observed_variance",
        &observation_intensity_observed_variance)
      .def("intensity_corrected_value",
        &observation_intensity_corrected_value)
      .def("intensity_corrected_variance",
        &observation_intensity_corrected_variance)
      .def("update_centroid_mm", 
        &observation_update_centroid_mm, (
          arg("detector"),
          arg("scan")));
  }

}}} // namespace dials::af::boost_python

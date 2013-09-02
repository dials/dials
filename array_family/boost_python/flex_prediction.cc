/*
 * flex_prediction.cc
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
#include <cctbx/miller.h>
#include <scitbx/array_family/boost_python/flex_wrapper.h>
#include <scitbx/array_family/ref_reductions.h>
#include <scitbx/array_family/boost_python/ref_pickle_double_buffered.h>
#include <scitbx/array_family/boost_python/flex_pickle_double_buffered.h>
#include <dials/model/data/prediction.h>
#include <dials/error.h>

namespace dials { namespace af { namespace boost_python {

  using namespace boost::python;
  using namespace scitbx::af::boost_python;
  using scitbx::af::const_ref;
  using scitbx::af::shared;
  using scitbx::vec2;
  using scitbx::vec3;
  using dials::model::Prediction;
  
  typedef cctbx::miller::index<> MillerIndex;

  static
  shared<MillerIndex> prediction_miller_index(
      const const_ref<Prediction> &obj) {
     shared<MillerIndex> result(obj.size());
     for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = obj[i].miller_index;
     }
     return result;
  }

  static
  shared< vec3<double> > prediction_beam_vector(
      const const_ref<Prediction> &obj) {
     shared< vec3<double> > result(obj.size());
     for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = obj[i].beam_vector;
     }
     return result;
  }

  static
  shared<int> prediction_panel(
      const const_ref<Prediction> &obj) {
     shared<int> result(obj.size());
     for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = obj[i].panel;
     }
     return result;
  }

  static
  shared<bool> prediction_entering(
      const const_ref<Prediction> &obj) {
     shared<bool> result(obj.size());
     for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = obj[i].entering;
     }
     return result;
  }

  static
  shared< vec3<double> > prediction_position_px(
      const const_ref<Prediction> &obj) {
     shared< vec3<double> > result(obj.size());
     for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = obj[i].position.px;
     }
     return result;
  }

  static
  shared< vec2<double> > prediction_position_px_coord(
      const const_ref<Prediction> &obj) {
     shared< vec2<double> > result(obj.size());
     for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = vec2<double>(obj[i].position.px[0], obj[i].position.px[1]);
     }
     return result;
  }

  static
  shared<double> prediction_position_frame(
      const const_ref<Prediction> &obj) {
     shared<double> result(obj.size());
     for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = obj[i].position.px[2];
     }
     return result;
  }

  static
  shared< vec3<double> > prediction_position_mm(
      const const_ref<Prediction> &obj) {
     shared< vec3<double> > result(obj.size());
     for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = obj[i].position.mm;
     }
     return result;
  }

  static
  shared< vec2<double> > prediction_position_mm_coord(
      const const_ref<Prediction> &obj) {
     shared< vec2<double> > result(obj.size());
     for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = vec2<double>(obj[i].position.mm[0], obj[i].position.mm[1]);
     }
     return result;
  }

  static
  shared<double> prediction_position_angle(
      const const_ref<Prediction> &obj) {
     shared<double> result(obj.size());
     for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = obj[i].position.mm[2];
     }
     return result;
  }

  struct prediction_to_string : pickle_double_buffered::to_string
  {
    using pickle_double_buffered::to_string::operator<<;

    prediction_to_string() {
      unsigned int version = 1;
      *this << version;
    }

    prediction_to_string& operator<<(const Prediction &val) {
      *this << val.miller_index[0]
            << val.miller_index[1]
            << val.miller_index[2]
            << val.beam_vector[0]
            << val.beam_vector[1]
            << val.beam_vector[2]
            << val.panel
            << val.entering
            << val.position.px[0]
            << val.position.px[1]
            << val.position.px[2]
            << val.position.mm[0]
            << val.position.mm[1]
            << val.position.mm[2];

      return *this;
    }
  };

  struct prediction_from_string : pickle_double_buffered::from_string
  {
    using pickle_double_buffered::from_string::operator>>;

    prediction_from_string(const char* str_ptr)
    : pickle_double_buffered::from_string(str_ptr) {
      *this >> version;
      DIALS_ASSERT(version == 1);
    }

    prediction_from_string& operator>>(Prediction &val) {
      *this >> val.miller_index[0]
            >> val.miller_index[1]
            >> val.miller_index[2]
            >> val.beam_vector[0]
            >> val.beam_vector[1]
            >> val.beam_vector[2]
            >> val.panel
            >> val.entering
            >> val.position.px[0]
            >> val.position.px[1]
            >> val.position.px[2]
            >> val.position.mm[0]
            >> val.position.mm[1]
            >> val.position.mm[2];

      return *this;
    }

    unsigned int version;
  };

  void export_flex_prediction()
  {
    scitbx::af::boost_python::flex_wrapper <
        Prediction, return_internal_reference<> >::plain("prediction")
      .def("miller_index",
        &prediction_miller_index)
      .def("beam_vector",
        &prediction_beam_vector)
      .def("panel",
        &prediction_panel)
      .def("entering",
        &prediction_entering)
      .def("position_px",
        &prediction_position_px)
      .def("position_px_coord",
        &prediction_position_px_coord)
      .def("position_frame",
        &prediction_position_frame)
      .def("position_mm",
        &prediction_position_mm)
      .def("position_mm_coord",
        &prediction_position_mm_coord)
      .def("position_angle",
        &prediction_position_angle)
      .def_pickle(flex_pickle_double_buffered<Prediction, 
        prediction_to_string, prediction_from_string>());
  }

}}} // namespace dials::af::boost_python

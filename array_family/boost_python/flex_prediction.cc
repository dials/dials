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
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/error.h>

namespace dials { namespace af { namespace boost_python {

  using namespace boost::python;
  using namespace scitbx::af::boost_python;
  using dials::model::Prediction;
  using scitbx::vec2;
  using scitbx::vec3;

  typedef cctbx::miller::index<> MillerIndex;

  /** @returns An array of miller indices */
  static af::shared<MillerIndex> prediction_miller_index(
    const af::const_ref<Prediction> &obj) {
    af::shared<MillerIndex> result(obj.size(), MillerIndex());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = obj[i].miller_index;
    }
    return result;
  }

  /** @returns An array of beam vectors */
  static af::shared<vec3<double> > prediction_beam_vector(
    const af::const_ref<Prediction> &obj) {
    af::shared<vec3<double> > result(obj.size(),
                                     af::init_functor_null<vec3<double> >());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = obj[i].beam_vector;
    }
    return result;
  }

  /** @returns An array of panel numbers */
  static af::shared<int> prediction_panel(const af::const_ref<Prediction> &obj) {
    af::shared<int> result(obj.size(), af::init_functor_null<int>());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = obj[i].panel;
    }
    return result;
  }

  /** @returns An array of entering flags */
  static af::shared<bool> prediction_entering(const af::const_ref<Prediction> &obj) {
    af::shared<bool> result(obj.size(), af::init_functor_null<bool>());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = obj[i].entering;
    }
    return result;
  }

  /** @returns An array of crystal ids */
  static af::shared<int> prediction_crystal(const af::const_ref<Prediction> &obj) {
    af::shared<int> result(obj.size(), af::init_functor_null<int>());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = obj[i].crystal;
    }
    return result;
  }

  /** @returns An array of pixel positions */
  static af::shared<vec3<double> > prediction_position_px(
    const af::const_ref<Prediction> &obj) {
    af::shared<vec3<double> > result(obj.size(),
                                     af::init_functor_null<vec3<double> >());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = obj[i].position.px;
    }
    return result;
  }

  /** @returns An array of pixel xy coordinates */
  static af::shared<vec2<double> > prediction_position_px_xy(
    const af::const_ref<Prediction> &obj) {
    af::shared<vec2<double> > result(obj.size(),
                                     af::init_functor_null<vec2<double> >());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = vec2<double>(obj[i].position.px[0], obj[i].position.px[1]);
    }
    return result;
  }

  /** @returns An array of frames */
  static af::shared<double> prediction_position_frame(
    const af::const_ref<Prediction> &obj) {
    af::shared<double> result(obj.size(), af::init_functor_null<double>());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = obj[i].position.px[2];
    }
    return result;
  }

  /** @returns An array of millimeter positions */
  static af::shared<vec3<double> > prediction_position_mm(
    const af::const_ref<Prediction> &obj) {
    af::shared<vec3<double> > result(obj.size(),
                                     af::init_functor_null<vec3<double> >());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = obj[i].position.mm;
    }
    return result;
  }

  /** @returns An array of millimeter xy coordinates */
  static af::shared<vec2<double> > prediction_position_mm_xy(
    const af::const_ref<Prediction> &obj) {
    af::shared<vec2<double> > result(obj.size(),
                                     af::init_functor_null<vec2<double> >());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = vec2<double>(obj[i].position.mm[0], obj[i].position.mm[1]);
    }
    return result;
  }

  /** @returns An array of angles */
  static af::shared<double> prediction_position_angle(
    const af::const_ref<Prediction> &obj) {
    af::shared<double> result(obj.size(), af::init_functor_null<double>());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = obj[i].position.mm[2];
    }
    return result;
  }

  /**
   * Class to convert predictions to strings for pickling
   */
  struct prediction_to_string : pickle_double_buffered::to_string {
    using pickle_double_buffered::to_string::operator<<;

    /** Initialise with the version */
    prediction_to_string() {
      unsigned int version = 1;
      *this << version;
    }

    /** Convert a single prediction instance to string */
    prediction_to_string &operator<<(const Prediction &val) {
      *this << val.miller_index[0] << val.miller_index[1] << val.miller_index[2]
            << val.beam_vector[0] << val.beam_vector[1] << val.beam_vector[2]
            << val.panel << val.entering << val.crystal << val.position.px[0]
            << val.position.px[1] << val.position.px[2] << val.position.mm[0]
            << val.position.mm[1] << val.position.mm[2];

      return *this;
    }
  };

  /**
   * Class to convert strings to predictions for unpickling
   */
  struct prediction_from_string : pickle_double_buffered::from_string {
    using pickle_double_buffered::from_string::operator>>;

    /** Initialise with the string */
    prediction_from_string(const char *str_ptr)
        : pickle_double_buffered::from_string(str_ptr) {
      *this >> version;
      DIALS_ASSERT(version == 1);
    }

    /** Convert a string to a single prediction instance */
    prediction_from_string &operator>>(Prediction &val) {
      *this >> val.miller_index[0] >> val.miller_index[1] >> val.miller_index[2]
        >> val.beam_vector[0] >> val.beam_vector[1] >> val.beam_vector[2] >> val.panel
        >> val.entering >> val.crystal >> val.position.px[0] >> val.position.px[1]
        >> val.position.px[2] >> val.position.mm[0] >> val.position.mm[1]
        >> val.position.mm[2];

      return *this;
    }

    unsigned int version;
  };

  void export_flex_prediction() {
    scitbx::af::boost_python::flex_wrapper<Prediction, return_internal_reference<> >::
      plain("prediction")
        .def("miller_index", &prediction_miller_index)
        .def("beam_vector", &prediction_beam_vector)
        .def("panel", &prediction_panel)
        .def("entering", &prediction_entering)
        .def("crystal", &prediction_crystal)
        .def("position_px", &prediction_position_px)
        .def("position_px_xy", &prediction_position_px_xy)
        .def("position_frame", &prediction_position_frame)
        .def("position_mm", &prediction_position_mm)
        .def("position_mm_xy", &prediction_position_mm_xy)
        .def("position_angle", &prediction_position_angle)
        .def_pickle(flex_pickle_double_buffered<Prediction,
                                                prediction_to_string,
                                                prediction_from_string>());
  }

}}}  // namespace dials::af::boost_python

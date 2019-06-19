/*
 * flex_intensity.cc
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
#include <scitbx/array_family/boost_python/ref_pickle_double_buffered.h>
#include <scitbx/array_family/boost_python/flex_pickle_double_buffered.h>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/model/data/observation.h>
#include <dials/error.h>

namespace dials { namespace af { namespace boost_python {

  using namespace boost::python;
  using namespace scitbx::af::boost_python;

  using dials::model::Intensity;
  using scitbx::vec2;
  using scitbx::vec3;

  /** @returns An array of observed intensity values */
  static af::shared<double> intensity_observed_value(
    const af::const_ref<Intensity> &obj) {
    af::shared<double> result(obj.size(), af::init_functor_null<double>());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = obj[i].observed.value;
    }
    return result;
  }

  /** @returns An array of observed intensity variances */
  static af::shared<double> intensity_observed_variance(
    const af::const_ref<Intensity> &obj) {
    af::shared<double> result(obj.size(), af::init_functor_null<double>());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = obj[i].observed.variance;
    }
    return result;
  }

  /** @returns An array of success values */
  static af::shared<bool> intensity_observed_success(
    const af::const_ref<Intensity> &obj) {
    af::shared<bool> result(obj.size(), af::init_functor_null<bool>());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = obj[i].observed.success;
    }
    return result;
  }

  /** @returns An array of corrected intensity values */
  static af::shared<double> intensity_corrected_value(
    const af::const_ref<Intensity> &obj) {
    af::shared<double> result(obj.size(), af::init_functor_null<double>());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = obj[i].corrected.value;
    }
    return result;
  }

  /** @returns An array of corrected intensity variances */
  static af::shared<double> intensity_corrected_variance(
    const af::const_ref<Intensity> &obj) {
    af::shared<double> result(obj.size(), af::init_functor_null<double>());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = obj[i].corrected.variance;
    }
    return result;
  }

  /** @returns An array of background intensity values */
  static af::shared<double> intensity_background_value(
    const af::const_ref<Intensity> &obj) {
    af::shared<double> result(obj.size(), af::init_functor_null<double>());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = obj[i].background.value;
    }
    return result;
  }

  /** @returns An array of background intensity variances */
  static af::shared<double> intensity_background_variance(
    const af::const_ref<Intensity> &obj) {
    af::shared<double> result(obj.size(), af::init_functor_null<double>());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = obj[i].background.variance;
    }
    return result;
  }

  void export_flex_intensity() {
    scitbx::af::boost_python::flex_wrapper<Intensity, return_internal_reference<> >::
      plain("intensity")
        .def("observed_value", &intensity_observed_value)
        .def("observed_variance", &intensity_observed_variance)
        .def("observed_success", &intensity_observed_success)
        .def("corrected_value", &intensity_corrected_value)
        .def("corrected_variance", &intensity_corrected_variance)
        .def("background_value", &intensity_background_value)
        .def("background_variance", &intensity_background_variance);
  }

}}}  // namespace dials::af::boost_python

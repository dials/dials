/*
 * flex_centroid.cc
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

  using dials::model::Centroid;
  using dxtbx::model::BeamBase;
  using dxtbx::model::Detector;
  using dxtbx::model::Scan;
  using scitbx::vec2;
  using scitbx::vec3;

  /** @returns An array of pixel positions */
  static af::shared<vec3<double> > centroid_px_position(
    const af::const_ref<Centroid> &obj) {
    af::shared<vec3<double> > result(obj.size(),
                                     af::init_functor_null<vec3<double> >());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = obj[i].px.position;
    }
    return result;
  }

  /** @returns An array of pixel variances */
  static af::shared<vec3<double> > centroid_px_variance(
    const af::const_ref<Centroid> &obj) {
    af::shared<vec3<double> > result(obj.size(),
                                     af::init_functor_null<vec3<double> >());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = obj[i].px.variance;
    }
    return result;
  }

  /** @returns An array of pixel standard error squareds */
  static af::shared<vec3<double> > centroid_px_std_err_sq(
    const af::const_ref<Centroid> &obj) {
    af::shared<vec3<double> > result(obj.size(),
                                     af::init_functor_null<vec3<double> >());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = obj[i].px.std_err_sq;
    }
    return result;
  }

  /** @returns An array of millimeter positions */
  static af::shared<vec3<double> > centroid_mm_position(
    const af::const_ref<Centroid> &obj) {
    af::shared<vec3<double> > result(obj.size(),
                                     af::init_functor_null<vec3<double> >());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = obj[i].mm.position;
    }
    return result;
  }

  /** @returns An array of millimeter variances */
  static af::shared<vec3<double> > centroid_mm_variance(
    const af::const_ref<Centroid> &obj) {
    af::shared<vec3<double> > result(obj.size(),
                                     af::init_functor_null<vec3<double> >());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = obj[i].mm.variance;
    }
    return result;
  }

  /** @returns An array of millimeter standard error squareds */
  static af::shared<vec3<double> > centroid_mm_std_err_sq(
    const af::const_ref<Centroid> &obj) {
    af::shared<vec3<double> > result(obj.size(),
                                     af::init_functor_null<vec3<double> >());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = obj[i].mm.std_err_sq;
    }
    return result;
  }

  /** @returns An array of pixel xy positions */
  static af::shared<vec2<double> > centroid_px_position_xy(
    const af::const_ref<Centroid> &obj) {
    af::shared<vec2<double> > result(obj.size(),
                                     af::init_functor_null<vec2<double> >());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = vec2<double>(obj[i].px.position[0], obj[i].px.position[1]);
    }
    return result;
  }

  /** @returns An array of frames */
  static af::shared<double> centroid_position_frame(
    const af::const_ref<Centroid> &obj) {
    af::shared<double> result(obj.size(), af::init_functor_null<double>());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = obj[i].px.position[2];
    }
    return result;
  }

  /** @returns An array of millimeter xy positions */
  static af::shared<vec2<double> > centroid_mm_position_xy(
    const af::const_ref<Centroid> &obj) {
    af::shared<vec2<double> > result(obj.size(),
                                     af::init_functor_null<vec2<double> >());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = vec2<double>(obj[i].mm.position[0], obj[i].mm.position[1]);
    }
    return result;
  }

  /** @returns An array of angles */
  static af::shared<double> centroid_position_angle(
    const af::const_ref<Centroid> &obj) {
    af::shared<double> result(obj.size(), af::init_functor_null<double>());
    for (std::size_t i = 0; i < result.size(); ++i) {
      result[i] = obj[i].mm.position[2];
    }
    return result;
  }

  /** Update the millimeter centroid positions of all observations */
  void centroid_update_mm(ref<Centroid> &obj, const Detector &d, const Scan &s) {
    for (std::size_t i = 0; i < obj.size(); ++i) {
      obj[i].update_mm(d, s);
    }
  }

  /** @returns The resolution of each observation */
  af::shared<double> centroid_resolution(af::ref<Centroid> &obj,
                                         std::size_t panel,
                                         const BeamBase &b,
                                         const Detector &d) {
    af::shared<double> result(obj.size(), af::init_functor_null<double>());
    for (std::size_t i = 0; i < obj.size(); ++i) {
      result[i] = obj[i].resolution(panel, b, d);
    }
    return result;
  }

  void export_flex_centroid() {
    scitbx::af::boost_python::flex_wrapper<Centroid, return_internal_reference<> >::
      plain("centroid")
        .def("px_position", &centroid_px_position)
        .def("px_variance", &centroid_px_variance)
        .def("px_std_err_eq", &centroid_px_std_err_sq)
        .def("mm_position", &centroid_mm_position)
        .def("mm_variance", &centroid_mm_variance)
        .def("mm_std_err_eq", &centroid_mm_std_err_sq)
        .def("px_position_xy", &centroid_px_position_xy)
        .def("position_frame", &centroid_position_frame)
        .def("mm_position_xy", &centroid_mm_position_xy)
        .def("position_angle", &centroid_position_angle)
        .def("update_mm",
             &centroid_update_mm,
             (boost::python::arg("detector"), boost::python::arg("scan")))
        .def("resolution",
             &centroid_resolution,
             (boost::python::arg("panel"),
              boost::python::arg("beam"),
              boost::python::arg("detector")));
  }

}}}  // namespace dials::af::boost_python

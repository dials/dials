/*
 * scan_varying_ray_predictor.cc
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
#include <dials/algorithms/spot_prediction/scan_varying_ray_predictor.h>
#include <dxtbx/model/scan.h>
#include <dxtbx/model/beam.h>
#include <dxtbx/model/goniometer.h>
#include <dxtbx/model/detector.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  static object predict1(const ScanVaryingRayPredictor &predictor,
                         const cctbx::miller::index<> &h,
                         const mat3<double> &A1,
                         const mat3<double> &A2,
                         int image,
                         std::size_t step) {
    boost::optional<Ray> result = predictor(h, A1, A2, image, step);
    return !result ? object() : object(*result);
  }

  static object predict2(const ScanVaryingRayPredictor &predictor,
                         const cctbx::miller::index<> &h,
                         const mat3<double> &A1,
                         const mat3<double> &A2,
                         const vec3<double> &s0a,
                         const vec3<double> &s0b,
                         int image,
                         std::size_t step) {
    boost::optional<Ray> result = predictor(h, A1, A2, s0a, s0b, image, step);
    return !result ? object() : object(*result);
  }

  void export_scan_varying_ray_predictor() {
    // Create and return the wrapper for the spot predictor object
    class_<ScanVaryingRayPredictor>("ScanVaryingRayPredictor", no_init)
      .def(init<vec3<double>, vec3<double>, int, vec2<double>, double>(
        (arg("s0"), arg("m2"), arg("frame0"), arg("dphi"), arg("dmin"))))
      .def("__call__",
           &predict1,
           (arg("hkl"), arg("A1"), arg("A2"), arg("image"), arg("step") = 1))
      .def("__call__",
           &predict2,
           (arg("hkl"),
            arg("A1"),
            arg("A2"),
            arg("s0a"),
            arg("s0b"),
            arg("image"),
            arg("step") = 1));
  }

}}}  // namespace dials::algorithms::boost_python

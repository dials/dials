/*
 * summation.cc
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
#include <dials/algorithms/integration/summation.h>

using namespace boost::python;

namespace dials { namespace algorithms { namespace boost_python {

  void export_summation()
  {
    class_ <IntegrateBySummation> ("IntegrateBySummation", no_init)
      .def(init <const flex_double&> ((
          arg("pixels"))))
      .def(init <const flex_double&, const flex_int&> ((
          arg("pixels"),
          arg("mask"))))
      .def(init <const flex_double&, const flex_vec3_double&> ((
          arg("pixels"),
          arg("points"))))
      .def("intensity", 
        &IntegrateBySummation::intensity)
      .def("variance",
        &IntegrateBySummation::variance)
      .def("standard_deviation",
        &IntegrateBySummation::standard_deviation)
      .def("centroid", 
        &IntegrateBySummation::centroid)
      .def("centroid_variance", 
        &IntegrateBySummation::centroid_variance)
      .def("centroid_standard_error_sq", 
        &IntegrateBySummation::centroid_standard_error_sq)
      .def("centroid_standard_error", 
        &IntegrateBySummation::centroid_standard_error)
      .def("centroid_covariance_matrix", 
        &IntegrateBySummation::centroid_covariance_matrix);
        
    Summation3d::integrator (Summation3d::*call_w_pixels)(
        const flex_int &pixels, 
        const flex_int &background,
        const flex_int &mask) const = 
      &Summation3d::operator();
    void (Summation3d::*call_w_reflection)(Reflection &r) const = 
      &Summation3d::operator();
    void (Summation3d::*call_w_reflection_list)(
        ReflectionList &reflections) const = 
      &Summation3d::operator();         
        
    class_<Summation3d>("Summation3dAlgorithm")
      .def("__call__", call_w_pixels)
      .def("__call__", call_w_reflection)
      .def("__call__", call_w_reflection_list);
  }

}}} // namespace = dials::algorithms::boost_python

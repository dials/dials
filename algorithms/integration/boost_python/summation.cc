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
      .def("intensity", &IntegrateBySummation::intensity)
      .def("centroid", &IntegrateBySummation::centroid)
      .def("variance", &IntegrateBySummation::variance)
      .def("standard_error_sq", &IntegrateBySummation::standard_error_sq)
      .def("standard_error", &IntegrateBySummation::standard_error)
      .def("covariance_matrix", &IntegrateBySummation::covariance_matrix);
  }

}}} // namespace = dials::algorithms::boost_python

/*
 * reflection_predictor.cc
 *
 *  Copyright (C) 2014 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <dials/algorithms/spot_prediction/reflection_predictor.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  void export_reflection_predictor()
  {
    class_<ScanStaticReflectionPredictor>("ScanStaticReflectionPredictor", no_init)
      .def("all_observable", &ScanStaticReflectionPredictor::all_observable)
      .def("observed", &ScanStaticReflectionPredictor::observed);
  }

}}} // namespace = dials::algorithms::boost_python

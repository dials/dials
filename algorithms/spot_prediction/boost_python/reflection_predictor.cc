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

    class_<ScanVaryingReflectionPredictor>("ScanVaryingReflectionPredictor", no_init)
      .def("all_observable", &ScanVaryingReflectionPredictor::all_observable);

    class_<StillsReflectionPredictor>("StillsReflectionPredictor", no_init)
      .def("observed", &StillsReflectionPredictor::observed);
  }

}}} // namespace = dials::algorithms::boost_python

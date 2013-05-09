/*
 * subtractor_strategy.cc
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
#include <dials/algorithms/background/subtractor_strategy.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  void export_subtractor_strategy()
  {
    void (SubtractorStrategy::*call_single)(Reflection&) const =
      &SubtractorStrategy::operator();
    flex_bool (SubtractorStrategy::*call_array)(ReflectionList&) const =
      &SubtractorStrategy::operator();
  
    class_<SubtractorStrategy, boost::noncopyable>(
        "SubtractorStrategy", no_init)
      .def("__call__", call_single)
      .def("__call__", call_array);
  }

}}} // namespace = dials::algorithms::boost_python

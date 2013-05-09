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
    class_<SubtractorStrategy, boost::noncopyable>(
        "SubtractorStrategy", no_init)
      .def("__call__", &SubtractorStrategy::operator());
  }

}}} // namespace = dials::algorithms::boost_python

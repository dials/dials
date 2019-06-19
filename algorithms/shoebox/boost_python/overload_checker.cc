/*
 * overload_checker.cc
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
#include <boost/python/iterator.hpp>
#include <dials/algorithms/shoebox/overload_checker.h>

namespace dials { namespace algorithms { namespace shoebox { namespace boost_python {

  using namespace boost::python;

  void export_overload_checker() {
    class_<OverloadChecker>("OverloadChecker")
      .def("add", &OverloadChecker::add)
      .def("__call__", &OverloadChecker::operator());
  }

}}}}  // namespace dials::algorithms::shoebox::boost_python

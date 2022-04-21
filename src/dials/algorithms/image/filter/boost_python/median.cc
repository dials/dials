/*
 * median.cc
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
#include <dials/algorithms/image/filter/median.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  template <typename T>
  void median_filter_suite() {
    def("median_filter", &median_filter<T>, (arg("image"), arg("kernel")));

    def("median_filter",
        &median_filter_masked<T>,
        (arg("image"), arg("mask"), arg("kernel"), arg("periodic") = false));
  }

  void export_median() {
    median_filter_suite<int>();
    median_filter_suite<float>();
    median_filter_suite<double>();
  }

}}}  // namespace dials::algorithms::boost_python

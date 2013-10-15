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

  void export_median() 
  {
    def("median_filter", &median_filter<int>, (
      arg("image"),
      arg("kernel")));
    def("median_filter", &median_filter<double>, (
      arg("image"),
      arg("kernel")));
    def("median_filter_masked", &median_filter_masked<int>, (
      arg("image"),
      arg("mask"),
      arg("kernel")));
    def("median_filter_masked", &median_filter_masked<double>, (
      arg("image"),
      arg("mask"),
      arg("kernel")));
  }

}}} // namespace = dials::algorithms::boost_python

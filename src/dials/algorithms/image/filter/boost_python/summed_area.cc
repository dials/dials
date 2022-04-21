/*
 * summed_area.cc
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
#include <dials/algorithms/image/filter/summed_area.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  template <typename T>
  void summed_area_suite() {
    def("summed_area_table", &summed_area_table<T>, (arg("image")));

    def("summed_area", &summed_area<T>, (arg("image"), arg("size")));
  }

  void export_summed_area() {
    summed_area_suite<int>();
    summed_area_suite<float>();
    summed_area_suite<double>();
  }

}}}  // namespace dials::algorithms::boost_python

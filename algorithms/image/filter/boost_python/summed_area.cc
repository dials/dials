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

  void export_summed_area() 
  {
    def("summed_area_table", &summed_area_table<int>, (arg("image")));
    def("summed_area_table", &summed_area_table<double>, (arg("image")));

    def("summed_area", &summed_area<int>, (arg("image"), arg("size")));
    def("summed_area", &summed_area<double>, (arg("image"), arg("size")));
  }

}}} // namespace = dials::algorithms::boost_python

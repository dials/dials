/*
 * mean_sdev_filter.cc
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
#include <boost_adaptbx/std_pair_conversion.h>
#include <dials/algorithms/peak_finding/mean_sdev_filter.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  void export_mean_sdev_filter() 
  {
    boost_adaptbx::std_pair_conversions::
      to_and_from_tuple<flex_double, flex_double>();  
  
    def("summed_area_table", &summed_area_table, (arg("image")));
    def("mean_filter", &mean_filter, (arg("image"), arg("size")));
    def("sdev_filter", &sdev_filter, (arg("image"), arg("mean"), arg("size")));
    def("mean_sdev_filter", &mean_sdev_filter, (arg("image"), arg("size")));
  }

}}} // namespace = dials::algorithms::boost_python

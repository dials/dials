/*
 * compute_centroid.cc
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
#include <dials/algorithms/centroid/compute_centroid.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  void export_compute_centroid()
  {
    void (ComputeCentroid::*call_array)(ReflectionList&) const = 
      &ComputeCentroid::operator();    
    
    void (ComputeCentroid::*call_single)(Reflection&) const = 
      &ComputeCentroid::operator();    

    class_<ComputeCentroid>("ComputeCentroid")
      .def("__call__", call_single)
      .def("__call__", call_array);
  }

}}} // namespace = dials::algorithms::boost_python

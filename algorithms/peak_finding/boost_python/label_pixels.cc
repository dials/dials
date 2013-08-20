/*
 * label_pixels.cc
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
#include <scitbx/array_family/boost_python/flex_wrapper.h>
#include <dials/algorithms/peak_finding/label_pixels.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  void export_label_pixels() 
  { 
    class_<LabelPixels>("LabelPixels", no_init)
      .def(init< vec3<int> >((arg("grid_size"))))
      .def("__call__", &LabelPixels::operator(), (arg("pixels")));
  }

}}} // namespace = dials::algorithms::boost_python

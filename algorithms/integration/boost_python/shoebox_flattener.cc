/*
 * shoebox_flattener.cc
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
#include <dials/algorithms/integration/shoebox_flattener.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  void export_shoebox_flattener()
  {
    class_<ShoeboxFlattener>("ShoeboxFlattener", no_init)
      .def(init<
          const GridSampler2D&,
          const af::const_ref< vec3<double> >&,
          const af::const_ref< Shoebox<> >&>())
      .def("index", &ShoeboxFlattener::index)
      .def("bbox", &ShoeboxFlattener::bbox)
      .def("data", &ShoeboxFlattener::data)
      .def("background", &ShoeboxFlattener::background)
      .def("mask", &ShoeboxFlattener::mask);
  }

}}} // namespace = dials::algorithms::boost_python


/*
 * shoebox.cc
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
#include <dials/algorithms/integration/shoebox.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  void export_shoebox()
  {
    class_ <Shoebox3d> ("Shoebox3d")
      .def(init<int6>())
      .def(init<int, int, int, int, int, int>())
      .def("min_x", &Shoebox3d::min_x)
      .def("max_x", &Shoebox3d::max_x)
      .def("min_y", &Shoebox3d::min_y)
      .def("max_y", &Shoebox3d::max_y)
      .def("min_z", &Shoebox3d::min_z)
      .def("max_z", &Shoebox3d::max_z)
      .def("range_x", &Shoebox3d::range_x)     
      .def("range_y", &Shoebox3d::range_y)
      .def("range_z", &Shoebox3d::range_z);
  }

}}} // namespace = dials::algorithms::boost_python

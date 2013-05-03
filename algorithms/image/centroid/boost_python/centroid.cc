/*
 * centroid2d.cc
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
#include <dials/algorithms/image/centroid/centroid.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  void export_centroid() 
  {
    class_<Centroid2d>("Centroid2d", no_init)
      .def(init<const flex_int&>((
        arg("image"))))
      .def("counts", &Centroid2d::counts)
      .def("position", &Centroid2d::position)
      .def("sq_width", &Centroid2d::sq_width)
      .def("variance", &Centroid2d::variance);

    class_<Centroid3d>("Centroid3d", no_init)
      .def(init<const flex_int&>((
        arg("image"))))
      .def("counts", &Centroid3d::counts)
      .def("position", &Centroid3d::position)
      .def("sq_width", &Centroid3d::sq_width)
      .def("variance", &Centroid3d::variance);
  }

}}} // namespace = dials::algorithms::boost_python

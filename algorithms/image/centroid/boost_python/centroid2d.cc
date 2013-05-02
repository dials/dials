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
#include <dials/algorithms/image/centroid/centroid2d.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  void export_centroid2d() 
  {
    class_<Centroid2d>("Centroid2d", no_init)
      .def(init<const flex_int&>((
        arg("image"))))
      .def("counts", &Centroid2d::counts)
      .def("position", &Centroid2d::position)
      .def("variance", &Centroid2d::variance)
      .def("variance_per_count", &Centroid2d::variance_per_count);

    class_<MaskedCentroid2d>("MaskedCentroid2d", no_init)
      .def(init<const flex_int&, const flex_int&>((
        arg("image"), arg("mask"))))
      .def("counts", &MaskedCentroid2d::counts)
      .def("position", &MaskedCentroid2d::position)
      .def("variance", &MaskedCentroid2d::variance)
      .def("variance_per_count", &MaskedCentroid2d::variance_per_count);
  }

}}} // namespace = dials::algorithms::boost_python

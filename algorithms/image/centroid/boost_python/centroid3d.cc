/*
 * centroid3d.cc
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
#include <dials/algorithms/image/centroid/centroid3d.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  void export_centroid3d() 
  {
    class_<Centroid3d>("Centroid3d", no_init)
      .def(init<const flex_int&>((
        arg("image"))))
      .def("counts", &Centroid3d::counts)
      .def("position", &Centroid3d::position)
      .def("variance", &Centroid3d::variance)
      .def("variance_per_count", &Centroid3d::variance_per_count);

    class_<MaskedCentroid3d>("MaskedCentroid3d", no_init)
      .def(init<const flex_int&, const flex_int&>((
        arg("image"), arg("mask"))))
      .def("counts", &MaskedCentroid3d::counts)
      .def("position", &MaskedCentroid3d::position)
      .def("variance", &MaskedCentroid3d::variance)
      .def("variance_per_count", &MaskedCentroid3d::variance_per_count);
  }

}}} // namespace = dials::algorithms::boost_python

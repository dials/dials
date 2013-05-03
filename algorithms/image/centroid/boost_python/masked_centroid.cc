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
#include <dials/algorithms/image/centroid/masked_centroid.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  void export_masked_centroid() 
  {
    class_<MaskedCentroid2d>("MaskedCentroid2d", no_init)
      .def(init<const flex_int&, const flex_int&>((
        arg("image"), arg("mask"))))
      .def("counts", &MaskedCentroid2d::counts)
      .def("position", &MaskedCentroid2d::position)
      .def("sq_width", &MaskedCentroid2d::sq_width)
      .def("variance", &MaskedCentroid2d::variance);  
  
    class_<MaskedCentroid3d>("MaskedCentroid3d", no_init)
      .def(init<const flex_int&, const flex_int&>((
        arg("image"), arg("mask"))))
      .def("counts", &MaskedCentroid3d::counts)
      .def("position", &MaskedCentroid3d::position)
      .def("sq_width", &MaskedCentroid3d::sq_width)
      .def("variance", &MaskedCentroid3d::variance);
  }

}}} // namespace = dials::algorithms::boost_python

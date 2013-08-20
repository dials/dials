/*
 * connected_components_ext.cc
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
#include <dials/algorithms/image/connected_components/connected_components.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  void export_connected_components()
  {
    class_<LabelImageStack>("LabelImageStack", no_init)
      .def(init<int2>((arg("size"))))
      .def("size", &LabelImageStack::size)
      .def("num_images", &LabelImageStack::num_images)
      .def("add_image", &LabelImageStack::add_image, (
        arg("image"), arg("mask")))
      .def("labels", &LabelImageStack::labels)
      .def("coords", &LabelImageStack::coords)
      .def("values", &LabelImageStack::values);
  }
  
  BOOST_PYTHON_MODULE(dials_algorithms_image_connected_components_ext)
  {
    export_connected_components();
  }

}}} // namespace = dials::algorithms::boost_python

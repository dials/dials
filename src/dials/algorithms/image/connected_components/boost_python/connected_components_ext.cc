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

  template <std::size_t DIM>
  void label_image_stack_wrapper(const char *name) {
    typedef LabelImageStack<DIM> label_type;
    class_<label_type>(name, no_init)
      .def(init<int2>((arg("size"))))
      .def("size", &label_type::size)
      .def("num_images", &label_type::num_images)
      .def("add_image", &label_type::add_image, (arg("image"), arg("mask")))
      .def("labels", &label_type::labels)
      .def("coords", &label_type::coords)
      .def("values", &label_type::values);
  }

  inline void label_pixels_wrapper(const char *name) {
    typedef LabelPixels label_type;
    class_<label_type>(name, no_init)
      .def(init<int3>((arg("size"))))
      .def("size", &label_type::size)
      .def("add_pixels", &label_type::add_pixels, (arg("values"), arg("coords")))
      .def("labels", &label_type::labels)
      .def("coords", &label_type::coords)
      .def("values", &label_type::values);
  }

  void export_connected_components() {
    label_image_stack_wrapper<2>("LabelImageStack2d");
    label_image_stack_wrapper<3>("LabelImageStack3d");
    label_pixels_wrapper("LabelPixels3d");
  }

  BOOST_PYTHON_MODULE(dials_algorithms_image_connected_components_ext) {
    export_connected_components();
  }

}}}  // namespace dials::algorithms::boost_python

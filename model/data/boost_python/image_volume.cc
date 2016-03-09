/*
 * image_volume.cc
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
#include <scitbx/array_family/flex_types.h>
#include <dials/model/data/image_volume.h>
#include <dials/config.h>

namespace dials { namespace model { namespace boost_python {

  using namespace boost::python;

  void export_image_volume()
  {
    class_<ImageVolume>("ImageVolume", no_init)
      .def(init<
          int,
          int,
          std::size_t,
          std::size_t>())
      .def("frame0", &ImageVolume::frame0)
      .def("frame1", &ImageVolume::frame1)
      .def("accessor", &ImageVolume::accessor)
      .def("data", &ImageVolume::data)
      .def("background", &ImageVolume::background)
      .def("mask", &ImageVolume::mask)
      .def("set_image", &ImageVolume::set_image<int>)
      .def("set_image", &ImageVolume::set_image<double>)
      ;

    class_<MultiPanelImageVolume>("MultiPanelImageVolume")
      .def("frame0", &MultiPanelImageVolume::frame0)
      .def("frame1", &MultiPanelImageVolume::frame1)
      .def("add", &MultiPanelImageVolume::add)
      .def("get", &MultiPanelImageVolume::get)
      .def("set_image", &MultiPanelImageVolume::set_image<int>)
      .def("set_image", &MultiPanelImageVolume::set_image<double>)
      .def("__len__", &MultiPanelImageVolume::size)
      ;
  }

}}} // namespace dials::model::boost_python


/*
 * pixel_list.cc
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
#include <dials/model/data/pixel_list.h>

namespace dials { namespace model { namespace boost_python {

  using namespace boost::python;

  void export_pixel_list()
  {
    class_<PixelList>("PixelList", no_init)
      .def(init<int2, int>((
        arg("size"),
        arg("start_frame"))))
      .def("size", &PixelList::size)
      .def("first_frame", &PixelList::first_frame)
      .def("last_frame", &PixelList::last_frame)
      .def("frame_range", &PixelList::frame_range)
      .def("num_frames", &PixelList::num_frames)
      .def("num_pixels", &PixelList::num_pixels)
      .def("add_image", &PixelList::add_image, (
        arg("image"), 
        arg("mask")))
      .def("coords", &PixelList::coords)
      .def("values", &PixelList::values)
      .def("labels_2d", &PixelList::labels_2d)
      .def("labels_3d", &PixelList::labels_3d);
  }

}}} // namespace dials::model::boost_python

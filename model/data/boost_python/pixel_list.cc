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

  struct PixelListPickleSuite : boost::python::pickle_suite {
    static boost::python::tuple getinitargs(const PixelList &obj) {
      return boost::python::make_tuple(
        obj.frame(), obj.size(), obj.value(), obj.index());
    }
  };

  void export_pixel_list() {
    class_<PixelList>("PixelList", no_init)
      .def(init<int,
                int2,
                const af::const_ref<double> &,
                const af::const_ref<std::size_t> &>(
        (arg("frame"), arg("size"), arg("value"), arg("index"))))
      .def(init<int,
                const af::const_ref<int, af::c_grid<2> > &,
                const af::const_ref<bool, af::c_grid<2> > &>(
        (arg("frame"), arg("size"), arg("value"), arg("index"))))
      .def(init<int,
                const af::const_ref<double, af::c_grid<2> > &,
                const af::const_ref<bool, af::c_grid<2> > &>(
        (arg("frame"), arg("size"), arg("value"), arg("index"))))
      .def("size", &PixelList::size)
      .def("frame", &PixelList::frame)
      .def("index", &PixelList::index)
      .def("value", &PixelList::value)
      .def("__len__", &PixelList::num_points)
      .def_pickle(PixelListPickleSuite());

    class_<PixelListLabeller>("PixelListLabeller")
      .def("add", &PixelListLabeller::add)
      .def("size", &PixelListLabeller::size)
      .def("num_pixels", &PixelListLabeller::num_pixels)
      .def("first_frame", &PixelListLabeller::first_frame)
      .def("last_frame", &PixelListLabeller::last_frame)
      .def("frame_range", &PixelListLabeller::frame_range)
      .def("num_frames", &PixelListLabeller::num_frames)
      .def("coords", &PixelListLabeller::coords)
      .def("values", &PixelListLabeller::values)
      .def("labels_3d", &PixelListLabeller::labels_3d)
      .def("labels_2d", &PixelListLabeller::labels_2d);
  }

}}}  // namespace dials::model::boost_python

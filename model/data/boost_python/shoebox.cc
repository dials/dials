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
#include <dials/model/data/shoebox.h>

namespace dials { namespace model { namespace boost_python {

  using namespace boost::python;

  int6 get_bbox(const Shoebox &obj) {
    return obj.bbox;
  }
  
  void set_bbox(Shoebox &obj, const int6 &v) {
    obj.bbox = v;
  }

  void export_shoebox()
  {
    class_<Shoebox>("Shoebox")
      .def_readwrite("data", &Shoebox::data)
      .def_readwrite("mask", &Shoebox::mask)
      .add_property("bbox", get_bbox, set_bbox)
      .def("allocate", &Shoebox::allocate)
      .def("deallocate", &Shoebox::deallocate)
      .def("xoffset", &Shoebox::xoffset)
      .def("yoffset", &Shoebox::yoffset)
      .def("zoffset", &Shoebox::zoffset)
      .def("offset", &Shoebox::offset)
      .def("xsize", &Shoebox::xsize)
      .def("ysize", &Shoebox::ysize)
      .def("zsize", &Shoebox::zsize)
      .def("size", &Shoebox::size)
      .def("consistent", &Shoebox::consistent);
  }

}}} // namespace dials::model::boost_python

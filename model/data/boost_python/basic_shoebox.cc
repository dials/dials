/*
 * basic_shoebox.cc
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
#include <dials/model/data/basic_shoebox.h>

namespace dials { namespace model { namespace boost_python {

  using namespace boost::python;
  using scitbx::af::flex_int;

  /** Set the data array as a flex array */
  static
  void set_data(BasicShoebox &obj, flex_int data) {
    DIALS_ASSERT(data.accessor().all().size() == 3);
    obj.data = af::versa<int, af::c_grid<3> >(
      data.handle(), af::c_grid<3>(data.accessor()));
  }

  static
  class_<BasicShoebox> basic_shoebox_wrapper(const char *name)
  {
    return class_<BasicShoebox>(name)
      .def(init<const BasicShoebox&>())
      .def(init<std::size_t, const int6&>())
      .add_property("data",
        make_getter(&BasicShoebox::data,
          return_value_policy<return_by_value>()),
        &set_data)
      .add_property("bbox",
        make_getter(&BasicShoebox::bbox,
          return_value_policy<return_by_value>()),
        make_setter(&BasicShoebox::bbox,
          return_value_policy<return_by_value>()))
      .def_readwrite("panel", &BasicShoebox::panel)
      .def("allocate", &BasicShoebox::allocate)
      .def("deallocate", &BasicShoebox::deallocate)
      .def("xoffset", &BasicShoebox::xoffset)
      .def("yoffset", &BasicShoebox::yoffset)
      .def("zoffset", &BasicShoebox::zoffset)
      .def("offset", &BasicShoebox::offset)
      .def("xsize", &BasicShoebox::xsize)
      .def("ysize", &BasicShoebox::ysize)
      .def("zsize", &BasicShoebox::zsize)
      .def("size", &BasicShoebox::size)
      .def("is_consistent", &BasicShoebox::is_consistent)
      .def("__eq__", &BasicShoebox::operator==)
      .def("__ne__", &BasicShoebox::operator!=);
  }

  void export_basic_shoebox()
  {
    basic_shoebox_wrapper("BasicShoebox");
  }

}}} // namespace dials::model::boost_python

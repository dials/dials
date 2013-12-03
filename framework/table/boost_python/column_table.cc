/*
 * column_table.cc
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
#include <dials/framework/table/column_table.h>

namespace dials { namespace framework { namespace boost_python {

  using namespace boost::python;

//  static
//  void get_data(const column_table &table, const std::string key) {

//  }

  void export_column_table() {
    class_<column_table>("column_table")
      .def("erase", &column_table::erase)
      .def("empty", &column_table::empty)
      .def("clear", &column_table::clear)
      .def("nrows", &column_table::nrows)
      .def("ncols", &column_table::ncols)
      .def("__len__", &column_table::size);
  }

}}} // namespace dials::framework::boost_python

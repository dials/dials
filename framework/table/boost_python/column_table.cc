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
#include <boost/python/slice.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/mpl/for_each.hpp>
#include <string>
#include <iterator>
#include <iostream>
#include <sstream>
#include <scitbx/array_family/flex_types.h>
#include <scitbx/boost_python/slice.h>
#include <dials/framework/table/column_table.h>
#include <dials/array_family/scitbx_shared_and_versa.h>
#include <dials/error.h>
#include <dials/framework/table/boost_python/column_table_suite.h>

namespace dials { namespace framework { namespace boost_python {

  using namespace boost::python;

  void export_column_table() {

    typedef column_type_generator<
      int,
      double,
      std::string
    >::type column_types;

    column_table_suite::column_data_wrapper<int>::wrap("column_data_int");
    column_table_suite::column_data_wrapper<double>::wrap("column_data_double");
    column_table_suite::column_data_wrapper<std::string>::wrap("column_data_std_string");
    column_table_suite::column_table_wrapper<column_types>::wrap("column_table");
  }

}}} // namespace dials::framework::boost_python

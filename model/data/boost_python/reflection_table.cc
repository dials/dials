/*
 * reflection_table.cc
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
#include <dials/framework/table/boost_python/column_table_suite.h>
#include <scitbx/array_family/tiny_types.h>
#include <scitbx/array_family/ref_reductions.h>
#include <scitbx/vec3.h>
#include <scitbx/vec2.h>

namespace dials { namespace model { namespace boost_python {

  using namespace boost::python;
  using namespace dials::framework;
  using namespace dials::framework::boost_python;
  using scitbx::vec2;
  using scitbx::vec3;
  using scitbx::af::int6;

  void export_reflection_table() {

    typedef column_type_generator<
      bool,
      int,
      std::size_t,
      double,
      std::string,
      vec2<double>,
      vec3<double>
    >::type column_types;

    column_table_suite::column_data_wrapper<bool>::wrap("column_bool");
    column_table_suite::column_data_wrapper<int>::wrap("column_int");
    column_table_suite::column_data_wrapper<std::size_t>::wrap("column_size_t");
    column_table_suite::column_data_wrapper<double>::wrap("column_double");
    column_table_suite::column_data_wrapper<std::string>::wrap("column_std_string");
    column_table_suite::column_data_wrapper< vec2<double> >::wrap("column_vec2_double");
    column_table_suite::column_data_wrapper< vec3<double> >::wrap("column_vec3_double");
    column_table_suite::column_table_wrapper<column_types>::wrap("reflection_table");
  }

}}} // namespace dials::model::boost_python

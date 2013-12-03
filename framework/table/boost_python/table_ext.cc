/*
 * table_ext.cc
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

namespace dials { namespace framework { namespace boost_python {

  using namespace boost::python;

  void export_column_table();

  BOOST_PYTHON_MODULE(dials_framework_table_ext)
  {
    export_column_table();
  }

}}} // namespace = dials::framework::boost_python

/*
 * data_ext.cc
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

namespace dials { namespace model { namespace boost_python {

  using namespace boost::python;

  void export_reflection();

  BOOST_PYTHON_MODULE(dials_model_data_ext)
  {
    export_reflection();
  }

}}} // namespace = dials::model::boost_python

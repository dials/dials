/*
 * preprocessor.cc
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
#include <dials/algorithms/integration/preprocessor.h>

using namespace boost::python;

namespace dials { namespace algorithms { namespace boost_python {

  void export_preprocessor() {

    class_<Preprocessor>("Preprocessor", no_init)
      .def(init<af::reflection_table>())
      .def("summary", &Preprocessor::summary)
      ;
  }

}}} // namespace = dials::algorithms::boost_python


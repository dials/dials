/*
 * helpers.cc
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
#include <dials/algorithms/background/helpers.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  void export_helpers()
  {
    def("set_shoebox_background_value",
      &set_shoebox_background_value, (
        arg("reflections"), arg("value")));
  }

}}} // namespace = dials::algorithms::boost_python

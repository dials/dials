/*
 * distance.cc
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
#include <dials/algorithms/image/filter/distance.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  void export_distance()
  {
    def("manhatten_distance",
      &manhatten_distance, (
        arg("data")));
  }

}}} // namespace = dials::algorithms::boost_python

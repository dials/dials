/*
 * unimodal.cc
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
#include <dials/algorithms/image/threshold/unimodal.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  void export_unimodal() {
    def("maximum_deviation", &maximum_deviation, (arg("histo")));
    def("probability_distribution",
        &probability_distribution,
        (arg("image"), arg("range")));
  }

}}}  // namespace dials::algorithms::boost_python

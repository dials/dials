/*
 * spot_finding_ext.cc
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
#include <dials/algorithms/spot_finding/helpers.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  BOOST_PYTHON_MODULE(dials_algorithms_spot_finding_ext) {
    class_<StrongSpotCombiner>("StrongSpotCombiner")
      .def("add", &StrongSpotCombiner::add)
      .def("shoeboxes", &StrongSpotCombiner::shoeboxes);
  }

}}}  // namespace dials::algorithms::boost_python

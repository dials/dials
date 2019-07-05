/*
 * scan_varying_helpers.cc
 *
 *  Copyright (C) 2014 STFC Rutherford Appleton Laboratory, UK.
 *
 *  Author: David Waterman
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <dials/algorithms/spot_prediction/scan_varying_helpers.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  void export_solve_quad() {
    def("solve_quad", reeke_detail::solve_quad, (arg("a"), arg("b"), arg("c")));
  }

}}}  // namespace dials::algorithms::boost_python

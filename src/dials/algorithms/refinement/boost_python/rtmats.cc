/*
 * rtmats.cc
 *
 *  Copyright (C) (2016) STFC Rutherford Appleton Laboratory, UK.
 *
 *  Author: David Waterman.
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include "../rtmats.h"

using namespace boost::python;

namespace dials { namespace refinement { namespace boost_python {

  void export_rtmats() {
    def("dR_from_axis_and_angle",
        &dR_from_axis_and_angle,
        (arg("axis"), arg("angle"), arg("deg") = false));
  }

}}}  // namespace dials::refinement::boost_python

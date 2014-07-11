/*
 * interpolate_profile2d.cc
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
#include <boost/python/iterator.hpp>
#include <dials/algorithms/integration/interpolate_profile2d.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  void export_interpolate_profile2d()
  {
    def("interpolate_profile2d", &interpolate_profile2d);
  }

}}} // namespace = dials::algorithms::boost_python

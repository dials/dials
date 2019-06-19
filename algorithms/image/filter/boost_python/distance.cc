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

  af::versa<int, af::c_grid<2> > manhattan_distance_wrapper(
    const af::const_ref<bool, af::c_grid<2> > &src,
    bool value) {
    af::versa<int, af::c_grid<2> > dst(src.accessor());
    manhattan_distance(src, value, dst.ref());
    return dst;
  }

  af::versa<int, af::c_grid<2> > chebyshev_distance_wrapper(
    const af::const_ref<bool, af::c_grid<2> > &src,
    bool value) {
    af::versa<int, af::c_grid<2> > dst(src.accessor());
    chebyshev_distance(src, value, dst.ref());
    return dst;
  }

  void export_distance() {
    def("manhattan_distance", &manhattan_distance_wrapper, (arg("data"), arg("value")));
    def("chebyshev_distance", &chebyshev_distance_wrapper, (arg("data"), arg("value")));
  }

}}}  // namespace dials::algorithms::boost_python

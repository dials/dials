/*
 * convolve.cc
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
#include <dials/algorithms/image/filter/convolve.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  template <typename FloatType>
  void convolve_suite() {
    def("convolve", &convolve<FloatType>, (arg("image"), arg("kernel")));

    def("convolve_row", &convolve_row<FloatType>, (arg("image"), arg("kernel")));

    def("convolve_col", &convolve_col<FloatType>, (arg("image"), arg("kernel")));
  }

  void export_convolve() {
    convolve_suite<float>();
    convolve_suite<double>();
  }

}}}  // namespace dials::algorithms::boost_python

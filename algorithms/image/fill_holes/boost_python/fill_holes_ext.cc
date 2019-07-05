/*
 * filter_ext.cc
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
#include <dials/algorithms/image/fill_holes/simple.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  BOOST_PYTHON_MODULE(dials_algorithms_image_fill_holes_ext) {
    def("simple_fill", &simple_fill, (arg("data"), arg("mask")));

    def(
      "diffusion_fill", &diffusion_fill, (arg("data"), arg("mask"), arg("niter") = 10));
  }

}}}  // namespace dials::algorithms::boost_python

/*
 * ext.cc
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
#include <dials/algorithms/background/median/creator.h>

namespace dials { namespace algorithms { namespace background { namespace boost_python {

  using namespace boost::python;

  BOOST_PYTHON_MODULE(dials_algorithms_background_median_ext) {
    def("create", &create_from_shoebox);
    def("create", &create_from_image_volume);
  }

}}}}  // namespace dials::algorithms::background::boost_python

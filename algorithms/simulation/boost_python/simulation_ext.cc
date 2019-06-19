/*
 * simulation_ext.cc
 *
 *  Copyright (C) 2014 Diamond Light Source
 *
 *  Author: James Parkhurst
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <dials/algorithms/simulation/reciprocal_space_helpers.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  BOOST_PYTHON_MODULE(dials_algorithms_simulation_ext) {
    def("simulate_reciprocal_space_gaussian", &simulate_reciprocal_space_gaussian);
    def("integrate_reciprocal_space_gaussian", &integrate_reciprocal_space_gaussian);
  }

}}}  // namespace dials::algorithms::boost_python

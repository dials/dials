/*
 * ideal_profile.cc
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
#include <dials/algorithms/reflection_basis/ideal_profile.h>

namespace dials { namespace algorithms { namespace reflection_basis {
  namespace transform { namespace boost_python {

  using namespace boost::python;

  void export_ideal_profile() {
    def("ideal_profile_float", &ideal_profile<float>);
    def("ideal_profile_double", &ideal_profile<double>);
  }

}}}}} // namespace dials::algorithms::reflexion_basis::transform::boost_python

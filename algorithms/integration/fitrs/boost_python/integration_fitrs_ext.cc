/*
 * integration_ext.cc
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

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  void export_profile_fitting_reciprocal_space();

  BOOST_PYTHON_MODULE(dials_algorithms_integration_fitrs_ext)
  {
    export_profile_fitting_reciprocal_space();
  }

}}} // namespace = dials::algorithms::boost_python

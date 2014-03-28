/*
 * profile_ext.cc
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
#include <iostream>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  void export_grid_sampler();
  void export_xds_circle_sampler();
  void export_reference_locator();
  void export_reference_learner();
  void export_fitting();

  BOOST_PYTHON_MODULE(dials_algorithms_integration_profile_ext)
  {
    export_grid_sampler();
    export_xds_circle_sampler();
    export_reference_locator();
    export_reference_learner();
    export_fitting();
  }

}}} // namespace = dials::algorithms::boost_python

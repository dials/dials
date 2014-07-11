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

  void export_luiso_s_2d_integration();
  void export_summation();
  void export_profile_fitting_reciprocal_space();
  void export_integrator_2d();
  void export_flattened_2d_profile_fitting();
  void export_interpolate_profile2d();

  BOOST_PYTHON_MODULE(dials_algorithms_integration_ext)
  {
    export_luiso_s_2d_integration();
    export_summation();
    export_profile_fitting_reciprocal_space();
    export_integrator_2d();
    export_flattened_2d_profile_fitting();
    export_interpolate_profile2d();
  }

}}} // namespace = dials::algorithms::boost_python

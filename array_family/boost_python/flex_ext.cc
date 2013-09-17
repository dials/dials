/*
 * flex_ext.cc
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
#include <scitbx/array_family/boost_python/c_grid_flex_conversions.h>
#include <dials/array_family/scitbx_shared_and_versa.h>

namespace dials { namespace af { namespace boost_python {

  using namespace boost::python;

  void export_flex_shoebox();
  void export_flex_centroid();
  void export_flex_intensity();
  void export_flex_observation();
  void export_flex_prediction();

  BOOST_PYTHON_MODULE(dials_array_family_flex_ext)
  {
    export_flex_shoebox();
    export_flex_centroid();
    export_flex_intensity();
    export_flex_observation();
    export_flex_prediction();
    
    scitbx::af::boost_python::c_grid_flex_conversions<double, af::c_grid<4> >();
  }

}}} // namespace = dials::af::boost_python

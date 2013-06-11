/*
 * background_ext.cc
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

  void export_helpers();
  void export_normal_discriminator();
  void export_poisson_discriminator();  
  void export_xds_subtractor();
  void export_fable_subtractor();
  
  BOOST_PYTHON_MODULE(dials_algorithms_background_ext)
  {
    export_helpers();
    export_normal_discriminator();
    export_poisson_discriminator();
    export_xds_subtractor();
    export_fable_subtractor();
  }

}}} // namespace = dials::algorithms::boost_python

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
#include <iostream>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  void export_sampler();
  void export_modeller();

  BOOST_PYTHON_MODULE(dials_algorithms_profile_model_modeller_ext) {
    export_sampler();
    export_modeller();
  }

}}}  // namespace dials::algorithms::boost_python

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

namespace dials { namespace algorithms { namespace background { namespace boost_python {

  using namespace boost::python;

  void export_outlier_rejector();
  void export_modeller();
  void export_creator();

  BOOST_PYTHON_MODULE(dials_algorithms_background_simple_ext) {
    export_outlier_rejector();
    export_modeller();
    export_creator();
  }

}}}}  // namespace dials::algorithms::background::boost_python

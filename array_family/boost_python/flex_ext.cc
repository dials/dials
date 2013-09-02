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

namespace dials { namespace af { namespace boost_python {

  using namespace boost::python;

  void export_flex_shoebox();
  void export_flex_observation();
  void export_flex_prediction();

  BOOST_PYTHON_MODULE(dials_array_family_flex_ext)
  {
    export_flex_shoebox();
    export_flex_observation();
    export_flex_prediction();
  }

}}} // namespace = dials::af::boost_python

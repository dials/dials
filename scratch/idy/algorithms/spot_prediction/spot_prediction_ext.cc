/*
 * spot_prediction_ext.cc
 *
 *  Copyright (C) 2015 Diamond Light Source
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#include <boost/python.hpp>
#include <boost/python/def.hpp>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  void export_experimental_reflection_predictor();

  BOOST_PYTHON_MODULE(idy_algorithms_spot_prediction_ext)
  {
    export_experimental_reflection_predictor();
  }

}}} // namespace = dials::algorithms::boost_python

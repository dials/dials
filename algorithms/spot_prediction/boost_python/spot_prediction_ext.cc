/*
 * spot_prediction_ext.cc
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

  void export_index_generator();
  void export_rotation_angles();
  void export_spot_predictor();
  void export_ray_predictor();

  BOOST_PYTHON_MODULE(dials_algorithms_spot_prediction_ext)
  {
    export_index_generator();
    export_rotation_angles();
    export_spot_predictor();
    export_ray_predictor();
  }

}}} // namespace = dials::algorithms::boost_python

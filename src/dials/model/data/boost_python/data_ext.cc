/*
 * data_ext.cc
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

namespace dials { namespace model { namespace boost_python {

  using namespace boost::python;

  void export_image_volume();
  void export_shoebox();
  void export_observation();
  void export_prediction();
  void export_pixel_list();
  void export_ray();
  void export_image();
  void export_adjacency_list();

  BOOST_PYTHON_MODULE(dials_model_data_ext) {
    export_image_volume();
    export_observation();
    export_prediction();
    export_shoebox();
    export_pixel_list();
    export_ray();
    export_image();
    export_adjacency_list();
  }

}}}  // namespace dials::model::boost_python

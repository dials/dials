/*
 * centroid_ext.cc
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

  void export_centroid_points();
  void export_centroid_image();
  void export_centroid_masked_image();

  BOOST_PYTHON_MODULE(dials_algorithms_image_centroid_ext)
  {
    export_centroid_points();
    export_centroid_image();
    export_centroid_masked_image();
  }

}}} // namespace = dials::algorithms::boost_python

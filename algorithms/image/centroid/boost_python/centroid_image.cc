/*
 * centroid_image.cc
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
#include <dials/algorithms/image/centroid/centroid_image.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  void export_centroid_image() 
  {
    class_<CentroidImage2d, bases<CentroidPoints <vec2<double> > > >(
        "CentroidImage2d", no_init)
      .def(init<const flex_double&>((arg("image"))));

    class_<CentroidImage3d, bases<CentroidPoints <vec3<double> > > >(
        "CentroidImage3d", no_init)
      .def(init<const flex_double&>((arg("image"))));
  }

}}} // namespace = dials::algorithms::boost_python

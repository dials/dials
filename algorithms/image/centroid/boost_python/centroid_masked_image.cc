/*
 * centroid_masked_image.cc
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
#include <dials/algorithms/image/centroid/centroid_masked_image.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  void export_centroid_masked_image() 
  {
    class_<CentroidMaskedImage2d, bases<CentroidPoints <vec2<double> > > >(
        "CentroidMaskedImage2d", no_init)
      .def(init<const af::const_ref<double, af::c_grid<2> >&, 
                const af::const_ref<int, af::c_grid<2> >&>((
        arg("image"), arg("mask"))));

    class_<CentroidMaskedImage3d, bases<CentroidPoints <vec3<double> > > >(
        "CentroidMaskedImage3d", no_init)
      .def(init<const af::const_ref<double, af::c_grid<3> >&, 
                const af::const_ref<int, af::c_grid<3> >&>((
        arg("image"), arg("mask"))));
  }

}}} // namespace = dials::algorithms::boost_python

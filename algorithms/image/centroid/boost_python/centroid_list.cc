/*
 * centroid_list.cc
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
#include <dials/algorithms/image/centroid/centroid_list.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  void export_centroid_list() 
  {
    class_<CentroidList2d>("CentroidList2d", no_init)
      .def(init<const flex_int&, const flex_vec2_double&>((
        arg("value"), arg("coord"))))
      .def("counts", &CentroidList2d::counts)
      .def("position", &CentroidList2d::position)
      .def("variance", &CentroidList2d::variance)
      .def("variance_per_count", &CentroidList2d::variance_per_count);

    class_<CentroidList3d>("CentroidList3d", no_init)
      .def(init<const flex_int&, const flex_vec3_double&>((
        arg("value"), arg("coord"))))
      .def("counts", &CentroidList3d::counts)
      .def("position", &CentroidList3d::position)
      .def("variance", &CentroidList3d::variance)
      .def("variance_per_count", &CentroidList3d::variance_per_count);
  }

}}} // namespace = dials::algorithms::boost_python

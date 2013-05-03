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
#include <scitbx/vec2.h>
#include <scitbx/vec3.h>
#include <dials/algorithms/image/centroid/centroid_list.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;
  using scitbx::vec2;
  using scitbx::vec3;

  template <typename CoordType>
  class_<CentroidList<CoordType> > class_centroid_list(const char *name) {
    typedef CentroidList<CoordType> CentroidListType;  
    
    return class_<CentroidListType>(name, no_init)
      .def(init<const flex_int&, const typename CentroidListType::flex_type&>((
        arg("pixels"), arg("coord"))))
      .def("counts", &CentroidListType::counts)
      .def("position", &CentroidListType::position)
      .def("sq_width", &CentroidListType::sq_width)
      .def("variance", &CentroidListType::variance);
  }
  
  CentroidList<vec2<double> > centroid_list_2d(const flex_int &pixels, 
      const flex<vec2<double> >::type &coords) {
    return CentroidList<vec2<double> >(pixels, coords);
  }

  CentroidList<vec3<double> > centroid_list_3d(const flex_int &pixels, 
      const flex<vec3<double> >::type &coords) {
    return CentroidList<vec3<double> >(pixels, coords);
  }

  void export_centroid_list() 
  {
    class_centroid_list<vec2<double> >("CentroidList2d");
    class_centroid_list<vec3<double> >("CentroidList3d");
    
    def("centroid_coord_list", &centroid_list_2d, (
      arg("pixels"), arg("coords")));
    def("centroid_coord_list", &centroid_list_3d, (
      arg("pixels"), arg("coords")));
  }

}}} // namespace = dials::algorithms::boost_python

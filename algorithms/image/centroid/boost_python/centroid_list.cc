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
      .def(init<const flex_double&, 
                const typename CentroidListType::flex_type&>((
          arg("pixels"), arg("coord"))))
      .def("sum_pixels", 
        &CentroidListType::sum_pixels)
      .def("sum_pixels_sq", 
        &CentroidListType::sum_pixels_sq)
      .def("sum_pixels_coords", 
        &CentroidListType::sum_pixels_coords)      
      .def("sum_pixels_delta_sq", 
        &CentroidListType::sum_pixels_delta_sq)
      .def("sum_pixels_delta_cross", 
        &CentroidListType::sum_pixels_delta_cross)
      .def("mean", 
        &CentroidListType::mean)
      .def("biased_variance", 
        &CentroidListType::biased_variance)
      .def("unbiased_variance",
        &CentroidListType::unbiased_variance)
      .def("biased_standard_error_sq", 
        &CentroidListType::biased_standard_error_sq)
      .def("unbiased_standard_error_sq", 
        &CentroidListType::unbiased_standard_error_sq)
      .def("covariance_matrix", 
        &CentroidListType::covariance_matrix);
  }
  
  CentroidList<vec2<double> > centroid_list_2d(const flex_double &pixels, 
      const flex<vec2<double> >::type &coords) {
    return CentroidList<vec2<double> >(pixels, coords);
  }

  CentroidList<vec3<double> > centroid_list_3d(const flex_double &pixels, 
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

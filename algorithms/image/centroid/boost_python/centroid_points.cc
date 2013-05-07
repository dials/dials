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
#include <dials/algorithms/image/centroid/centroid_points.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;
  using scitbx::vec2;
  using scitbx::vec3;

  template <typename CoordType>
  class_<CentroidPoints<CoordType> > class_centroid_points(const char *name) {
    typedef CentroidPoints<CoordType> CentroidPointsType;  
    
    return class_<CentroidPointsType>(name, no_init)
      .def(init<const flex_double&, 
                const typename CentroidPointsType::flex_type&>((
          arg("pixels"), arg("coord"))))
      .def("sum_pixels", 
        &CentroidPointsType::sum_pixels)
      .def("sum_pixels_sq", 
        &CentroidPointsType::sum_pixels_sq)
      .def("sum_pixels_coords", 
        &CentroidPointsType::sum_pixels_coords)      
      .def("sum_pixels_delta_sq", 
        &CentroidPointsType::sum_pixels_delta_sq)
      .def("sum_pixels_delta_cross", 
        &CentroidPointsType::sum_pixels_delta_cross)
      .def("mean", 
        &CentroidPointsType::mean)
      .def("biased_variance", 
        &CentroidPointsType::biased_variance)
      .def("unbiased_variance",
        &CentroidPointsType::unbiased_variance)
      .def("biased_standard_error_sq", 
        &CentroidPointsType::biased_standard_error_sq)
      .def("unbiased_standard_error_sq", 
        &CentroidPointsType::unbiased_standard_error_sq)
      .def("covariance_matrix", 
        &CentroidPointsType::covariance_matrix);
  }
  
  CentroidPoints<vec2<double> > centroid_points_2d(const flex_double &pixels, 
      const flex<vec2<double> >::type &coords) {
    return CentroidPoints<vec2<double> >(pixels, coords);
  }

  CentroidPoints<vec3<double> > centroid_points_3d(const flex_double &pixels, 
      const flex<vec3<double> >::type &coords) {
    return CentroidPoints<vec3<double> >(pixels, coords);
  }

  void export_centroid_points() 
  {
    class_centroid_points<vec2<double> >("CentroidPoints2d");
    class_centroid_points<vec3<double> >("CentroidPoints3d");
    
    def("centroid_points", &centroid_points_2d, (
      arg("pixels"), arg("coords")));
    def("centroid_points", &centroid_points_3d, (
      arg("pixels"), arg("coords")));
  }

}}} // namespace = dials::algorithms::boost_python

/*
 * multi_plane_geometry.cc
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
#include <boost/format.hpp>
#include <string>
#include <dials/model/experiment/multi_plane_geometry.h>

namespace dials { namespace model { namespace boost_python {

  using namespace boost::python;

  void export_beam_multi_plane_intersection()
  {
    typedef BeamMultiPlaneIntersection::coord_type coord_type;
    
    coord_type (BeamMultiPlaneIntersection::*call_w_beam_vector)(
      vec3 <double>) const = &BeamMultiPlaneIntersection::operator();
    vec2 <double> (BeamMultiPlaneIntersection::*call_w_beam_vector_and_plane)(
      vec3 <double>, int) const = &BeamMultiPlaneIntersection::operator();
  
    // Export the object to calculate beam->plane intersections
    class_<BeamMultiPlaneIntersection>("BeamMultiPlaneIntersection", no_init)
      .def(init<const flex_mat3_double&, const flex_double4 &>((
        arg("D"), arg("extents"))))
      .def("__call__", 
        call_w_beam_vector, (
          arg("s1")))
      .def("__call__", 
        call_w_beam_vector_and_plane, (
          arg("s1"), arg("plane")));    
  }
  
  void export_multi_plane_to_lab_transform()
  {
    typedef MultiPlaneToLabTransform::coord_type coord_type;
    
    vec3 <double> (MultiPlaneToLabTransform::*call_w_plane_xy)(
      int, vec2 <double>) const = &MultiPlaneToLabTransform::operator();
    vec3 <double> (MultiPlaneToLabTransform::*call_w_coord_type)(
      coord_type) const = &MultiPlaneToLabTransform::operator();
      
    // Export the object to calculate the lab from pane coordinates
    class_<MultiPlaneToLabTransform>("MultiPlaneToLabTransform", no_init)
      .def(init<const flex_mat3_double &>((
        arg("d"))))
      .def("__call__",
        call_w_plane_xy, (
          arg("plane"), arg("xy")))
      .def("__call__",
        call_w_coord_type, (
          arg("pxy")));
  }

  void export_multi_plane_geometry()
  {
    export_beam_multi_plane_intersection();
    export_multi_plane_to_lab_transform();
  }  

}}} // namespace = dials::model::boost_python

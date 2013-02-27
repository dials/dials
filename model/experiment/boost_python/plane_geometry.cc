/*
 * plane_geometry.cc
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
#include <dials/model/experiment/plane_geometry.h>

namespace dials { namespace model { namespace boost_python {

  using namespace boost::python;

  void export_plane_geometry()
  {
    // Export the object to calculate beam->plane intersections
    class_<BeamPlaneIntersection>("BeamPlaneIntersection", no_init)
      .def(init<mat3<double> >((
        arg("D"))))
      .def("__call__", 
        &BeamPlaneIntersection::operator(), (
          arg("s1")));
          
    // Export the object to calculate the lab from pane coordinates
    class_<PlaneToLabTransform>("PlaneToLabTransform", no_init)
      .def(init<mat3<double> >((
        arg("d"))))
      .def("__call__",
        &PlaneToLabTransform::operator(), (
          arg("xy")));
  }
  

}}} // namespace = dials::model::boost_python

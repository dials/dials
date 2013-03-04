/*
 * ray_intersector.cc
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
#include <dxtbx/model/detector.h>
#include <dials/model/data/reflection.h>
#include <dials/algorithms/spot_prediction/ray_intersector.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  void export_ray_intersector()
  {
    // Typedef the different overloads for operator()
    Reflection (RayIntersector::*intersect_single)(const Reflection&) const = 
      &RayIntersector::operator();
    ReflectionList (RayIntersector::*intersect_array)(
      const ReflectionList&) const = &RayIntersector::operator();

    // Typedef the different overloads for operator()
    Reflection (RayIntersector::*intersect_single_w_panel)(const Reflection&, 
      int panel) const = &RayIntersector::operator();
    ReflectionList (RayIntersector::*intersect_array_w_panel)(
      const ReflectionList&, int panel) const = &RayIntersector::operator();

    // Create and return the wrapper for the spot predictor object
    class_ <RayIntersector> ("RayIntersector", no_init)
      .def(init <const Detector &> ((
        arg("detector"))))
      .def("__call__", intersect_single, (
        arg("reflection")))
      .def("__call__", intersect_array, (
        arg("reflections")))
      .def("__call__", intersect_single_w_panel, (
        arg("reflection")))
      .def("__call__", intersect_array_w_panel, (
        arg("reflections")));
  }

}}} // namespace = dials::spot_prediction::boost_python

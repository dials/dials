/*
 * ray_intersection.cc
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
#include <dials/algorithms/spot_prediction/ray_intersection.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  void export_ray_intersection()
  {
    // Pointers to function overloads
    Reflection (*ray_intersection_single)(const Detector&, const Reflection&)  
      = &ray_intersection;
    shared<Reflection> (*ray_intersection_array)(const Detector&, 
      const ReflectionList&) = &ray_intersection;
    Reflection (*ray_intersection_single_w_panel)(const Detector&, 
      const Reflection&, std::size_t) = &ray_intersection;
    shared<Reflection> (*ray_intersection_array_w_panel)(const Detector&, 
      const ReflectionList&, std::size_t) = &ray_intersection;

    // Export all the ray intersection functions
    def("ray_intersection", ray_intersection_single, 
      (arg("detector"), arg("reflection")));
    def("ray_intersection", ray_intersection_single_w_panel, 
      (arg("detector"), arg("reflection"), arg("panel")));
    def("ray_intersection", ray_intersection_array, 
      (arg("detector"), arg("reflection_list")));
    def("ray_intersection", ray_intersection_array_w_panel, 
      (arg("detector"), arg("reflection_list"), arg("panel")));
  }

}}} // namespace = dials::spot_prediction::boost_python

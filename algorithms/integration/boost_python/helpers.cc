/*
 * helpers.cc
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
#include <dials/algorithms/integration/helpers.h>

namespace dials { namespace algorithms { namespace boost_python { 

  using namespace boost::python;
      
  void export_helpers() 
  {
    void (*filter_reflection_by_detector_mask)(Reflection&, const flex_bool&,
        int2) = &filter_by_detector_mask;
  
    void (*filter_reflection_list_by_detector_mask)(ReflectionList&, 
      const flex_bool&, int2) = &filter_by_detector_mask;
  
    def("is_bbox_outside_image_range", &is_bbox_outside_image_range);
    def("does_bbox_contain_bad_pixels", &does_bbox_contain_bad_pixels);
    def("is_bbox_valid", &is_bbox_valid);
    def("filter_by_detector_mask", filter_reflection_by_detector_mask);
    def("filter_by_detector_mask", filter_reflection_list_by_detector_mask);
    def("filter_by_bbox_volume", filter_by_bbox_volume);
  }

}}} // namespace = dials::algorithms::boost_python


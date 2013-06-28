/*
 * reflection_profile_helpers.cc
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
#include <dials/algorithms/integration/reflection_profile_helpers.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  void export_reflection_profile_helpers()
  {
    def("allocate_reflection_profiles",
      &allocate_reflection_profiles, (
        arg("reflections"),
        arg("shoebox_default") = 0,
        arg("shoebox_mask_default") = 1,
        arg("shoebox_background_default") = 0));
  
    def("copy_single_image_pixels",
      &copy_single_image_pixels, (
        arg("image"), 
        arg("array_index"), 
        arg("reflection_indices"), 
        arg("reflections"),
        arg("gain_map"),
        arg("dark_map")));
        
    def("construct_image_mask_from_shoeboxes", 
      &construct_image_mask_from_shoeboxes, (
        arg("mask"), 
        arg("array_index"),
        arg("index"),
        arg("reflections"),
        arg("kernel_size"))); 
  }

}}} // namespace = dials::algorithms::boost_python

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
    def("assign_strong_spots",
      &assign_strong_spots, (
        arg("image"), 
        arg("array_index"), 
        arg("reflection_indices"), 
        arg("reflections")));        
  }

}}} // namespace = dials::algorithms::boost_python

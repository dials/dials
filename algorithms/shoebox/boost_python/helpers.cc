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
#include <dials/algorithms/shoebox/helpers.h>

namespace dials { namespace algorithms { namespace shoebox {
  namespace boost_python {

  using namespace boost::python;

  void export_helpers()
  {
    def("allocate", (void(*)(ReflectionList&))&allocate, (
      arg("reflection_list")));
    def("allocate", (void(*)(ReflectionList&,int))&allocate, (
      arg("reflection_list"), arg("mask_default")));
    def("deallocate", &deallocate, (
      arg("reflection_list")));
    def("assign_strong_spots",
      &assign_strong_spots, (
        arg("image"), 
        arg("array_index"), 
        arg("reflection_indices"), 
        arg("reflections")));        
  }

}}}} // namespace = dials::algorithms::shoebox::boost_python

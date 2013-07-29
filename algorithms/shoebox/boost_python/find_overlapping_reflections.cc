/*
 * find_overlapping_reflections.cc
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
#include <boost/python/iterator.hpp>
#include <boost_adaptbx/std_pair_conversion.h>
#include <dials/algorithms/shoebox/find_overlapping_reflections.h>

namespace dials { namespace algorithms { namespace shoebox {
  namespace boost_python {

  using namespace boost::python;

  void export_find_overlapping_reflections()
  {
    def("find_overlapping_reflections", 
      &find_overlapping_reflections, (
        arg("reflection_list")));
  }

}}}} // namespace = dials::algorithms::shoebox::boost_python

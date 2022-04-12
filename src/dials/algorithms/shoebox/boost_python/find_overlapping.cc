/*
 * find_overlapping.cc
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
#include <dials/algorithms/shoebox/find_overlapping.h>

namespace dials { namespace algorithms { namespace shoebox { namespace boost_python {

  using namespace boost::python;

  void export_find_overlapping() {
    def("find_overlapping", &find_overlapping, (arg("bboxes")));
    def("find_overlapping", &find_overlapping_multi_panel, (arg("bbox"), arg("panel")));

    class_<OverlapFinder>("OverlapFinder").def("__call__", &OverlapFinder::operator());
  }

}}}}  // namespace dials::algorithms::shoebox::boost_python

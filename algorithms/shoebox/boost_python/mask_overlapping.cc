/*
 * mask_overlapping.cc
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
#include <dials/algorithms/shoebox/mask_overlapping.h>

namespace dials { namespace algorithms { namespace shoebox { namespace boost_python {

  using namespace boost::python;

  void export_mask_overlapping() {
    class_<MaskOverlapping>("MaskOverlapping")
      .def("__call__",
           &MaskOverlapping::operator(),
           (arg("shoeboxes"), arg("coords"), arg("adjacency_list")));
  }

}}}}  // namespace dials::algorithms::shoebox::boost_python

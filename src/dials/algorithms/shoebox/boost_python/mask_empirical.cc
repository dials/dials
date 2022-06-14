/*
 * mask_empirical.cc
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
#include <dials/algorithms/shoebox/mask_empirical.h>

namespace dials { namespace algorithms { namespace shoebox { namespace boost_python {

  using namespace boost::python;

  void export_mask_empirical() {
    class_<MaskEmpirical>("MaskEmpirical", no_init)
      .def(init<const af::reflection_table&>((arg("reference"))))
      .def("__call__", &MaskEmpirical::mask);
  }

}}}}  // namespace dials::algorithms::shoebox::boost_python

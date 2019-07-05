/*
 * index_generator.cc
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
#include <dials/algorithms/spot_prediction/index_generator.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  void export_index_generator() {
    class_<IndexGenerator>("IndexGenerator")
      .def(init<cctbx::uctbx::unit_cell const&,
                cctbx::sgtbx::space_group_type const&,
                double>(
        (arg("unit_cell"), arg("space_group_type"), arg("resolution_d_min"))))
      .def("next", &IndexGenerator::next)
      .def("to_array", &IndexGenerator::to_array);
  }

}}}  // namespace dials::algorithms::boost_python

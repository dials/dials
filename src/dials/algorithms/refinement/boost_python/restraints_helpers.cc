/*
 * restraints_helpers.cc
 *
 *  Copyright (C) (2016) STFC Rutherford Appleton Laboratory, UK.
 *
 *  Author: David Waterman.
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include "../restraints/restraints_helpers.h"

using namespace boost::python;

namespace dials { namespace refinement { namespace boost_python {

  void export_calculate_cell_gradients() {
    class_<CalculateCellGradients>("CalculateCellGradients", no_init)
      .def(init<mat3<double>, af::const_ref<mat3<double> > >((arg("B"), arg("dB_dp"))))
      .def("da_dp", &CalculateCellGradients::da_dp)
      .def("db_dp", &CalculateCellGradients::db_dp)
      .def("dc_dp", &CalculateCellGradients::dc_dp)
      .def("daa_dp", &CalculateCellGradients::daa_dp)
      .def("dbb_dp", &CalculateCellGradients::dbb_dp)
      .def("dcc_dp", &CalculateCellGradients::dcc_dp);
  }

}}}  // namespace dials::refinement::boost_python

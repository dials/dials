/*
 * obs_manager_helpers.cc
 *
 *  Copyright (C) (2016) STFC Rutherford Appleton Laboratory, UK.
 *
 *  Author: David Waterman
 *
 *  This code is distributed under the BSD license, a copy of which is
 *  included in the root directory of this package.
 */
#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <dials/algorithms/scaling/scale_parameterisation_helpers.h>

namespace dials { namespace scaling { namespace boost_python {

  using namespace boost::python;

  void export_row_multiply()
  {
    def("row_multiply", &row_multiply, (
      arg("m"),
      arg("v")));
  }

}}} // namespace = dials::algorithms::boost_python


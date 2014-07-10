/*
 * flattened_2d_profile_fitting.cc
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
#include <dials/algorithms/integration/flattened_2d_profile_fitting.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  void export_flattened_2d_profile_fitting()
  {
    class_<Flattened2DProfileFitting>(
        "Flattened2DProfileFitting", no_init)
      .def("intensity", &Flattened2DProfileFitting::intensity)
      .def("variance", &Flattened2DProfileFitting::variance )
      ;
  }

}}} // namespace = dials::algorithms::boost_python


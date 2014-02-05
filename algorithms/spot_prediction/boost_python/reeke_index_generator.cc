/*
 * reeke_index_generator.cc
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
#include <dials/algorithms/spot_prediction/reeke_index_generator.h>

namespace dials { namespace algorithms { namespace boost_python {

  using namespace boost::python;

  void export_reeke_index_generator()
  {
    class_<ReekeIndexGenerator>("ReekeIndexGenerator", no_init)
      .def(init<mat3<double>,
                mat3<double>,
                vec3<double>,
                vec3<double>,
                double,
                int>())
      .def("next", &ReekeIndexGenerator::next)
      .def("to_array", &ReekeIndexGenerator::to_array);

  }

}}} // namespace = dials::algorithms::boost_python

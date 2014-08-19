/*
 * interface.cc
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
#include <dials/algorithms/integration/interface.h>

using namespace boost::python;

namespace dials { namespace algorithms { namespace boost_python {

  void export_interface() {

    class_<IntegrationManagerData3D>("IntegrationManagerData3D", no_init)
      .def(init<af::reflection_table,
                vec2<double>,
                vec2<int>,
                double>((
          arg("reflections"),
          arg("oscillation"),
          arg("array_range"),
          arg("block_size"))))
      .def("__len__", &IntegrationManagerData3D::size)
      .def("finished", &IntegrationManagerData3D::finished)
      .def("block", &IntegrationManagerData3D::block)
      .def("to_process", &IntegrationManagerData3D::to_process)
      .def("to_include", &IntegrationManagerData3D::to_include)
      .def("to_not_process", &IntegrationManagerData3D::to_not_process)
      .def("split", &IntegrationManagerData3D::split)
      .def("accumulate", &IntegrationManagerData3D::accumulate)
      .def("data", &IntegrationManagerData3D::data)
      ;
  }

}}} // namespace = dials::algorithms::boost_python

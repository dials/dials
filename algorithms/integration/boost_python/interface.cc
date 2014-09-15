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

    class_<IntegrationManagerExecutor>("IntegrationManagerExecutor", no_init)
      .def(init<af::reflection_table,
                vec2<int>,
                double>((
          arg("reflections"),
          arg("array_range"),
          arg("block_size"))))
      .def("__len__", &IntegrationManagerExecutor::size)
      .def("finished", &IntegrationManagerExecutor::finished)
      .def("accumulate", &IntegrationManagerExecutor::accumulate)
      .def("split", &IntegrationManagerExecutor::split)
      .def("job", &IntegrationManagerExecutor::job)
      .def("data", &IntegrationManagerExecutor::data)
      .def("ignored", &IntegrationManagerExecutor::ignored)
      ;
  }

}}} // namespace = dials::algorithms::boost_python

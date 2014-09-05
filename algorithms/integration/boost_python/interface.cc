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

    class_<IntegrationManager3DExecutor>("IntegrationManager3DExecutor", no_init)
      .def(init<af::reflection_table,
                vec2<int>,
                double>((
          arg("reflections"),
          arg("array_range"),
          arg("block_size"))))
      .def("__len__", &IntegrationManager3DExecutor::size)
      .def("finished", &IntegrationManager3DExecutor::finished)
      .def("accumulate", &IntegrationManager3DExecutor::accumulate)
      .def("split", &IntegrationManager3DExecutor::split)
      .def("job", &IntegrationManager3DExecutor::job)
      .def("data", &IntegrationManager3DExecutor::data)
      .def("ignored", &IntegrationManager3DExecutor::ignored)
      ;
  }

}}} // namespace = dials::algorithms::boost_python

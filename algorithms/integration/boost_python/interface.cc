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

    class_<IntegrationJobCalculator>("IntegrationJobCalculator", no_init)
      .def(init< vec2<int>, double >((
          arg("array_range"),
          arg("block_size"))))
      .def("jobs", &IntegrationJobCalculator::jobs)
      ;

    class_<IntegrationManagerExecutor>("IntegrationManagerExecutor", no_init)
      .def(init<const IntegrationJobCalculator&,
                af::reflection_table>((
          arg("jobcalculator"),
          arg("reflections"))))
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

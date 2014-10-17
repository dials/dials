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

    class_<JobCalculator>("JobCalculator", no_init)
      .def(init< vec2<int>, double >((
          arg("array_range"),
          arg("block_size"))))
      .def("jobs", &JobCalculator::jobs)
      ;

    class_<ReflectionManager>("ReflectionManager", no_init)
      .def(init<const JobCalculator&,
                af::reflection_table>((
          arg("jobcalculator"),
          arg("reflections"))))
      .def("__len__", &ReflectionManager::size)
      .def("finished", &ReflectionManager::finished)
      .def("accumulate", &ReflectionManager::accumulate)
      .def("split", &ReflectionManager::split)
      .def("job", &ReflectionManager::job)
      .def("data", &ReflectionManager::data)
      ;
  }

}}} // namespace = dials::algorithms::boost_python

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

    class_<GroupList::Group>("Group", no_init)
      .def("index", &GroupList::Group::index)
      .def("nindex", &GroupList::Group::nindex)
      .def("expr", &GroupList::Group::expr)
      .def("nexpr", &GroupList::Group::nexpr)
      .def("frames", &GroupList::Group::frames)
      .def("nframes", &GroupList::Group::nframes)
      ;

    class_<GroupList>("GroupList")
      .def("__len__", &GroupList::size)
      .def("__getitem__", &GroupList::operator[],
          return_internal_reference<>())
      ;

    class_<JobList::Job>("Job", no_init)
      .def("index", &JobList::Job::index)
      .def("expr", &JobList::Job::expr)
      .def("nexpr", &JobList::Job::nexpr)
      .def("frames", &JobList::Job::frames)
      .def("nframes", &JobList::Job::nframes)
      ;

    class_<JobList>("JobList")
      .def("add", &JobList::add)
      .def("__len__", &JobList::size)
      .def("__getitem__", &JobList::operator[],
          return_internal_reference<>())
      ;

    class_<ReflectionManager>("ReflectionManager", no_init)
      .def(init<const JobList&,
                af::reflection_table>((
          arg("jobs"),
          arg("data"))))
      .def("__len__", &ReflectionManager::size)
      .def("finished", &ReflectionManager::finished)
      .def("accumulate", &ReflectionManager::accumulate)
      .def("split", &ReflectionManager::split)
      .def("job", &ReflectionManager::job,
          return_internal_reference<>())
      .def("data", &ReflectionManager::data)
      ;
  }

}}} // namespace = dials::algorithms::boost_python

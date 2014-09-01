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

  struct callback_helper {
    object obj_;
    callback_helper(object obj) : obj_(obj) {}
    af::reflection_table operator()(af::reflection_table data) const {
      return extract<af::reflection_table>(obj_(data));
    }
  };

  IntegrationTask3DExecutor* make_integration_task_3d_executor(
      const IntegrationTask3DSpec &spec,
      object callback) {
    return new IntegrationTask3DExecutor(spec, callback_helper(callback));
  }

  void export_interface() {

    class_<IntegrationTask3DSpec>("IntegrationTask3DSpec", no_init)
      .def(init< af::reflection_table,
                 std::size_t,
                 const af::const_ref< tiny<int,2> >&,
                 const af::const_ref< std::size_t >&,
                 const af::const_ref< std::size_t >&,
                 const af::const_ref< bool >&>())
      .def("data", &IntegrationTask3DSpec::data)
      .def("job", &IntegrationTask3DSpec::job)
      .def("frame0", &IntegrationTask3DSpec::frame0)
      .def("frame1", &IntegrationTask3DSpec::frame1)
      .def("nframes", &IntegrationTask3DSpec::nframes)
      .def("njobs", &IntegrationTask3DSpec::njobs)
      ;

    class_<IntegrationTask3DExecutor>("IntegrationTask3DExecutor", no_init)
      .def("__init__", make_constructor(
            &make_integration_task_3d_executor))
      .def("next", &IntegrationTask3DExecutor::next)
      .def("frame0", &IntegrationTask3DExecutor::frame0)
      .def("frame1", &IntegrationTask3DExecutor::frame1)
      .def("frame", &IntegrationTask3DExecutor::frame)
      .def("nframes", &IntegrationTask3DExecutor::nframes)
      .def("finished", &IntegrationTask3DExecutor::finished)
      ;

    class_<IntegrationManager3DExecutor>("IntegrationManager3DExecutor", no_init)
      .def(init<af::reflection_table,
                vec2<int>,
                double,
                std::size_t>((
          arg("reflections"),
          arg("array_range"),
          arg("block_size"),
          arg("num_tasks"))))
      .def("__len__", &IntegrationManager3DExecutor::size)
      .def("finished", &IntegrationManager3DExecutor::finished)
      .def("task", &IntegrationManager3DExecutor::task)
      .def("split", &IntegrationManager3DExecutor::split)
      .def("accumulate", &IntegrationManager3DExecutor::accumulate)
      .def("data", &IntegrationManager3DExecutor::data)
      .def("jobs", &IntegrationManager3DExecutor::jobs);
      ;

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

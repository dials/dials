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
      af::reflection_table data,
      const af::const_ref< tiny<int,2> > &jobs,
      std::size_t npanels,
      object callback) {
    return new IntegrationTask3DExecutor(data, jobs, npanels, callback_helper(callback));
  }

  void export_interface() {

    class_<IntegrationTask3DExecutor>("IntegrationTask3DExecutor", no_init)
      .def("__init__", make_constructor(
            &make_integration_task_3d_executor))
      .def("next", &IntegrationTask3DExecutor::next)
      .def("frame0", &IntegrationTask3DExecutor::frame0)
      .def("frame1", &IntegrationTask3DExecutor::frame1)
      .def("frame", &IntegrationTask3DExecutor::frame)
      .def("nframes", &IntegrationTask3DExecutor::nframes)
      .def("finished", &IntegrationTask3DExecutor::finished)
      .def("data", &IntegrationTask3DExecutor::data)
      ;

    class_<IntegrationTask3DMultiExecutor>("IntegrationTask3DMultiExecutorBase", no_init)
      .def(init<
          af::reflection_table,
          tiny<int,2>,
          std::size_t>())
      .def("next", &IntegrationTask3DMultiExecutor::next)
      .def("frame0", &IntegrationTask3DMultiExecutor::frame0)
      .def("frame1", &IntegrationTask3DMultiExecutor::frame1)
      .def("frame", &IntegrationTask3DMultiExecutor::frame)
      .def("nframes", &IntegrationTask3DMultiExecutor::nframes)
      .def("finished", &IntegrationTask3DMultiExecutor::finished)
      .def("data", &IntegrationTask3DMultiExecutor::data)
      ;

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

    class_<IntegrationManager3DMultiExecutor>("IntegrationManager3DMultiExecutor", no_init)
      .def(init<af::reflection_table,
                vec2<int>,
                double>((
          arg("reflections"),
          arg("array_range"),
          arg("block_size"))))
      .def("__len__", &IntegrationManager3DMultiExecutor::size)
      .def("finished", &IntegrationManager3DMultiExecutor::finished)
      .def("accumulate", &IntegrationManager3DMultiExecutor::accumulate)
      .def("split", &IntegrationManager3DMultiExecutor::split)
      .def("job", &IntegrationManager3DMultiExecutor::job)
      .def("data", &IntegrationManager3DMultiExecutor::data)
      .def("ignored", &IntegrationManager3DMultiExecutor::ignored)
      ;
  }

}}} // namespace = dials::algorithms::boost_python

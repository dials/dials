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

    class_<Summary::ImageSummary>("ImageSummary", no_init)
      .def("full", &Summary::ImageSummary::full)
      .def("part", &Summary::ImageSummary::part)
      .def("sum_ios", &Summary::ImageSummary::sum_ios)
      .def("prf_ios", &Summary::ImageSummary::prf_ios)
      ;

    class_<Summary::ResolutionSummary>("ResolutionSummary", no_init)
      .def("nbins", &Summary::ResolutionSummary::nbins)
      .def("bin_range", &Summary::ResolutionSummary::bin_range)
      .def("bins", &Summary::ResolutionSummary::bins)
      .def("sum_ios", &Summary::ResolutionSummary::sum_ios)
      .def("sum_ios_full", &Summary::ResolutionSummary::sum_ios_full)
      .def("sum_ios_part", &Summary::ResolutionSummary::sum_ios_part)
      .def("prf_ios", &Summary::ResolutionSummary::prf_ios)
      .def("prf_ios_full", &Summary::ResolutionSummary::prf_ios_full)
      .def("prf_ios_part", &Summary::ResolutionSummary::prf_ios_part)
      ;

    class_<Summary::WholeSummary>("WholeSummary", no_init)
      .def("sum_ios", &Summary::WholeSummary::sum_ios)
      .def("sum_ios_full", &Summary::WholeSummary::sum_ios_full)
      .def("sum_ios_part", &Summary::WholeSummary::sum_ios_part)
      .def("prf_ios", &Summary::WholeSummary::prf_ios)
      .def("prf_ios_full", &Summary::WholeSummary::prf_ios_full)
      .def("prf_ios_part", &Summary::WholeSummary::prf_ios_part)
      ;

    class_<Summary>("Summary", no_init)
      .def(init<af::reflection_table,
                int2,
                std::size_t>())
      .def("image_summary", &Summary::image_summary)
      .def("resolution_summary", &Summary::resolution_summary)
      .def("whole_summary", &Summary::whole_summary)
      ;
  }

}}} // namespace = dials::algorithms::boost_python

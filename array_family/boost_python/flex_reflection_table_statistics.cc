/*
 * reflection_table_statistics.cc
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
#include <dials/array_family/reflection_table_statistics.h>

using namespace boost::python;

namespace dials { namespace af { namespace boost_python {

  void image_summary_wrapper() {

    typedef Summary::ImageSummary SummaryClass;

    class_<SummaryClass::Element>("ImageSummaryElement", no_init)
      .def_readonly("frame", &SummaryClass::Element::frame)
      .def_readonly("full", &SummaryClass::Element::full)
      .def_readonly("part", &SummaryClass::Element::part)
      .add_property("sum_ios", &SummaryClass::Element::sum_ios_value)
      .add_property("prf_ios", &SummaryClass::Element::prf_ios_value)
      ;

    class_<SummaryClass::Data>("ImageSummaryData", no_init)
      .def("__getitem__", &SummaryClass::Data::at)
      .def("__len__", &SummaryClass::Data::size)
      ;

    class_<SummaryClass>("ImageSummary", no_init)
      .def("data", &SummaryClass::data)
      .def("__len__", &SummaryClass::size)
      ;
  }

  void resolution_summary_wrapper() {

    typedef Summary::ResolutionSummary SummaryClass;

    class_<SummaryClass::Element>("ResolutionSummaryElement", no_init)
      .add_property("d", &SummaryClass::Element::d_value)
      .add_property("sum_ios", &SummaryClass::Element::sum_ios_value)
      .add_property("sum_ios_full", &SummaryClass::Element::sum_ios_full_value)
      .add_property("sum_ios_part", &SummaryClass::Element::sum_ios_part_value)
      .add_property("prf_ios", &SummaryClass::Element::prf_ios_value)
      .add_property("prf_ios_full", &SummaryClass::Element::prf_ios_full_value)
      .add_property("prf_ios_part", &SummaryClass::Element::prf_ios_part_value)
      ;

    class_<SummaryClass::Data>("ResolutionSummaryData", no_init)
      .def("__getitem__", &SummaryClass::Data::at)
      .def("__len__", &SummaryClass::Data::size)
      ;

    class_<SummaryClass>("ResolutionSummary", no_init)
      .def("data", &SummaryClass::data)
      .def("__len__", &SummaryClass::size)
      ;
  }

  void whole_summary_wrapper() {

    typedef Summary::WholeSummary SummaryClass;

    class_<SummaryClass::Data>("WholeSummaryData", no_init)
      .add_property("sum_ios", &SummaryClass::Data::sum_ios_value)
      .add_property("sum_ios_full", &SummaryClass::Data::sum_ios_full_value)
      .add_property("sum_ios_part", &SummaryClass::Data::sum_ios_part_value)
      .add_property("prf_ios", &SummaryClass::Data::prf_ios_value)
      .add_property("prf_ios_full", &SummaryClass::Data::prf_ios_full_value)
      .add_property("prf_ios_part", &SummaryClass::Data::prf_ios_part_value)
      ;

    class_<SummaryClass>("WholeSummary", no_init)
      .def("data", &SummaryClass::data)
      .def("__len__", &SummaryClass::size)
      ;
  }

  void summary_wrapper() {

    image_summary_wrapper();
    resolution_summary_wrapper();
    whole_summary_wrapper();

    class_<Summary>("Summary", no_init)
      .def(init<af::reflection_table,
                const af::const_ref<int2>&,
                std::size_t>())
      .def("image_summary", &Summary::image_summary)
      .def("resolution_summary", &Summary::resolution_summary)
      .def("whole_summary", &Summary::whole_summary)
      ;
  }

  void export_flex_reflection_table_statistics() {
    summary_wrapper();

/*     class_<Summary::ImageSummary>("ImageSummary", no_init) */
/*       .def("full", &Summary::ImageSummary::full) */
/*       .def("part", &Summary::ImageSummary::part) */
/*       .def("sum_ios", &Summary::ImageSummary::sum_ios) */
/*       .def("prf_ios", &Summary::ImageSummary::prf_ios) */
/*       ; */

/*     class_<Summary::ResolutionSummary>("ResolutionSummary", no_init) */
/*       .def("nbins", &Summary::ResolutionSummary::nbins) */
/*       .def("bin_range", &Summary::ResolutionSummary::bin_range) */
/*       .def("bins", &Summary::ResolutionSummary::bins) */
/*       .def("sum_ios", &Summary::ResolutionSummary::sum_ios) */
/*       .def("sum_ios_full", &Summary::ResolutionSummary::sum_ios_full) */
/*       .def("sum_ios_part", &Summary::ResolutionSummary::sum_ios_part) */
/*       .def("prf_ios", &Summary::ResolutionSummary::prf_ios) */
/*       .def("prf_ios_full", &Summary::ResolutionSummary::prf_ios_full) */
/*       .def("prf_ios_part", &Summary::ResolutionSummary::prf_ios_part) */
/*       ; */

/*     class_<Summary::WholeSummary>("WholeSummary", no_init) */
/*       .def("sum_ios", &Summary::WholeSummary::sum_ios) */
/*       .def("sum_ios_full", &Summary::WholeSummary::sum_ios_full) */
/*       .def("sum_ios_part", &Summary::WholeSummary::sum_ios_part) */
/*       .def("prf_ios", &Summary::WholeSummary::prf_ios) */
/*       .def("prf_ios_full", &Summary::WholeSummary::prf_ios_full) */
/*       .def("prf_ios_part", &Summary::WholeSummary::prf_ios_part) */
/*       ; */

/*     class_<Summary>("Summary", no_init) */
/*       .def(init<af::reflection_table, */
/*                 int2, */
/*                 std::size_t>()) */
/*       .def("image_summary", &Summary::image_summary) */
/*       .def("resolution_summary", &Summary::resolution_summary) */
/*       .def("whole_summary", &Summary::whole_summary) */
/*       ; */
    }

}}} // namespace dials::af::boost_python
